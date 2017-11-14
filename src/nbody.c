// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#include "nbody.h"
#include "verlet_ode_solver.h"
#include "barnes_hut_tree.h"

typedef struct
{
  // Parallel stuff.
  MPI_Comm comm;
  int rank, nprocs;

  // Gravitational coupling constant.
  real_t G;

  // Initial data.
  body_array_t* bodies;

  // Barnes-Hut algorithm stuff.
  barnes_hut_tree_t* tree;

  // Time integrator
  ode_solver_t* solver;

  // Solution data.
  real_t* U;
  real_t v_min, a_max;
  real_t E0; // initial energy.
} nbody_t;

static real_t nbody_total_energy(nbody_t* nb)
{
  // The total energy in the system is the sum of gravitational 
  // potentials plus kinetic energies.
  real_t E = 0.0, G = nb->G;
  for (size_t b = 0; b < nb->bodies->size; ++b)
  {
    body_t* b1 = nb->bodies->data[b];
    real_t m1 = b1->m;
    real_t vx = nb->U[6*b+3];
    real_t vy = nb->U[6*b+4];
    real_t vz = nb->U[6*b+5];
    E += 0.5 * m1 * (vx*vx + vy*vy + vz*vz);
    point_t x1 = {.x = nb->U[6*b], .y = nb->U[6*b+1], .z = nb->U[6*b+2]};
    for (size_t bb = 0; bb < nb->bodies->size; ++bb)
    {
      if (b == bb) continue;
      body_t* b2 = nb->bodies->data[bb];
      real_t m2 = b2->m;
      point_t x2 = {.x = nb->U[6*bb], .y = nb->U[6*bb+1], .z = nb->U[6*bb+2]};
      real_t r2 = point_square_distance(&x1, &x2);
      E += G * m1 * m2 / r2;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &E, 1, MPI_REAL_T, MPI_SUM, nb->comm);
  return E;
}

static void nbody_init(void* context, real_t t)
{
  nbody_t* nb = context;

  // Populate the solution vector.
  nb->v_min = REAL_MAX;
  for (size_t b = 0; b < nb->bodies->size; ++b)
  {
    body_t* body = nb->bodies->data[b];
    nb->U[6*b]   = body->x.x;
    nb->U[6*b+1] = body->x.y;
    nb->U[6*b+2] = body->x.z;
    nb->U[6*b+3] = body->v.x;
    nb->U[6*b+4] = body->v.y;
    nb->U[6*b+5] = body->v.z;
    nb->v_min = MIN(nb->v_min, vector_mag(&body->v));
  }

  // Compute the initial energy of the system.
  nb->E0 = nbody_total_energy(nb);
}

static real_t nbody_max_dt(void* context, real_t t, char* reason)
{
//  nbody_t* nb = context;

  // The timestep is limited by the ratio of the maximum acceleration to 
  // the minimum speed.
  snprintf(reason, POLYMEC_MODEL_MAXDT_REASON_SIZE,
           "ratio of max acceleration to minimum speed");
  return REAL_MAX;//0.2 * (nb->a_max / nb->v_min);
}

static real_t nbody_advance(void* context, real_t max_dt, real_t t)
{
  nbody_t* nb = context;
  polymec_suspend_fpe();
  bool solved = ode_solver_advance(nb->solver, t, t + max_dt, nb->U);
  polymec_restore_fpe();
  if (!solved)
    return 0.0;
  else
    return max_dt;
}

static void nbody_finalize(void* context, int step, real_t t)
{
  // Compute the final energy and compare it to the initial energy.
  nbody_t* nb = context;
  real_t E = nbody_total_energy(nb);
  log_detail("Fractional energy change: %g", (E - nb->E0) / nb->E0 - 1.0);
}

static void nbody_dtor(void* context)
{
  nbody_t* nb = context;
  body_array_free(nb->bodies);
  polymec_free(nb->U);
  ode_solver_free(nb->solver);
  if (nb->tree != NULL)
    barnes_hut_tree_free(nb->tree);
  polymec_free(nb);
}

// Solver functions.
static int brute_force_accel(void* context, real_t t, real_t* U, real_t* dvdt)
{
  nbody_t* nb = context;
  real_t G = nb->G;
  nb->a_max = 0.0;
  nb->v_min = REAL_MAX;

  // Extract point coordinates and masses from solution data.
  int N = (int)(nb->bodies->size);
  point_t x[N];
  real_t m[N];
  for (int i = 0; i < N; ++i)
  {
    x[i].x = U[6*i];
    x[i].y = U[6*i+1];
    x[i].z = U[6*i+2];
    m[i] = nb->bodies->data[i]->m;
  }

  // How many points on other processes?
  int Np[nb->nprocs];
  MPI_Allgather(&N, 1, MPI_INT, 
                Np, 1, MPI_INT, 
                nb->comm);
  int Ntot = 0;
  for (int p = 0; p < nb->nprocs; ++p)
    Ntot += Np[p];

  // Get point data from other processes.
  point_t* all_points = polymec_malloc(sizeof(point_t) * Ntot);
  MPI_Allgather(x, 3*Ntot, MPI_REAL_T, 
                all_points, 3*Ntot, MPI_REAL_T, 
                nb->comm);

  // Now compute forces on all pairs, with i ranging over the entirety of 
  // the points, and j over just the local ones.
  vector_t* all_accels = polymec_malloc(sizeof(vector_t) * Ntot);
  memset(all_accels, 0, sizeof(vector_t) * Ntot);
  for (size_t i = 0; i < Ntot; ++i)
  {
    point_t xi = all_points[i];
    vector_t a = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (size_t j = 0; j < N; ++j)
    {
      if (i == j) continue;
      real_t mj = m[j];
      point_t xj = x[j];
      vector_t rij;
      point_displacement(&xi, &xj, &rij);
      real_t r = vector_mag(&rij);
      real_t r3_inv = 1.0 / (r*r*r);
      a.x += G * mj * rij.x * r3_inv;
      a.y += G * mj * rij.y * r3_inv;
      a.z += G * mj * rij.z * r3_inv;
    }

    // Update the acceleration.
    all_accels[i] = a;
  }

  // Sum together all accelerations on all points over all processes.
  MPI_Allreduce(MPI_IN_PLACE, all_accels, 3*Ntot, 
                MPI_REAL_T, MPI_SUM, nb->comm);

  // Extract acceleration data for this process.
  int start = 0;
  for (int p = 0; p < nb->rank; ++p)
    start += Np[p];
  memcpy(dvdt, &(all_accels[start]), sizeof(vector_t) * N);

  // Clean up.
  polymec_free(all_accels);
  polymec_free(all_points);

  return 0;
}

static int barnes_hut_accel(void* context, real_t t, real_t* U, real_t* dvdt)
{
  nbody_t* nb = context;
  int N = (int)(nb->bodies->size);
  point_t x[N];
  real_t m[N];
  for (int i = 0; i < N; ++i)
  {
    x[i].x = U[6*i];
    x[i].y = U[6*i+1];
    x[i].z = U[6*i+2];
    m[i] = nb->bodies->data[i]->m;
  }
  barnes_hut_tree_compute_forces(nb->tree, nb->G, x, m, N, (vector_t*)dvdt);
  for (int i = 0; i < N; ++i)
  {
    dvdt[3*i] /= m[i];
    dvdt[3*i+1] /= m[i];
    dvdt[3*i+2] /= m[i];
  }

  return 0;
}

static body_array_t* partition_bodies(nbody_t* nb, body_array_t* bodies)
{
  START_FUNCTION_TIMER();

  int N = (int)(bodies->size);

  // Partition the bodies and their data. Since each process has all of the 
  // data, we just compute a partition vector.
  point_cloud_t* cloud = point_cloud_new(nb->comm, N);
  for (int b = 0; b < N; ++b)
    cloud->points[b] = bodies->data[b]->x;
  int64_t* P = partition_vector_from_point_cloud(cloud, nb->comm, NULL, 1.05);

  // Use the partition vector to cull the off-process bodies.
  body_array_t* local_bodies = body_array_new();
  for (int i = 0; i < N; ++i)
  {
    if (P[i] == nb->rank)
      body_array_append(local_bodies, body_clone(bodies->data[i]));
  }
  nb->bodies = local_bodies;

  // Clean up.
  polymec_free(P);
  point_cloud_free(cloud);

  STOP_FUNCTION_TIMER();
  return local_bodies;
}

//------------------------------------------------------------------------

model_t* brute_force_nbody_new(real_t G,
                               body_array_t* bodies)
{
  ASSERT(G >= 0.0);

  nbody_t* nb = polymec_malloc(sizeof(nbody_t));
  nb->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(nb->comm, &(nb->rank));
  MPI_Comm_size(nb->comm, &(nb->nprocs));
  nb->G = G;
  nb->bodies = partition_bodies(nb, bodies);
  nb->tree = NULL;
  nb->U = polymec_malloc(sizeof(real_t) * 6 * bodies->size);

  // Dispose of the original global list of bodies.
  body_array_free(bodies);

  // Set up a Verlet solver.
  nb->solver = verlet_ode_solver_new(nb->comm, (int)nb->bodies->size,
                                     nb, brute_force_accel, NULL);

  // Now create the model.
  model_vtable vtable = {.init = nbody_init,
                         .max_dt = nbody_max_dt,
                         .advance = nbody_advance,
                         .finalize = nbody_finalize,
                         .dtor = nbody_dtor};
  return model_new("Brute-Force N-body", nb, vtable, MODEL_SERIAL);
}

model_t* barnes_hut_nbody_new(real_t G,
                              real_t theta,
                              body_array_t* bodies)
{
  ASSERT(G >= 0.0);
  ASSERT(theta >= 0.0);

  nbody_t* nb = polymec_malloc(sizeof(nbody_t));
  nb->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(nb->comm, &(nb->rank));
  MPI_Comm_size(nb->comm, &(nb->nprocs));
  nb->G = G;
  nb->bodies = partition_bodies(nb, bodies);
  nb->tree = barnes_hut_tree_new(nb->comm, theta);
  nb->U = polymec_malloc(sizeof(real_t) * 6 * bodies->size);

  // Dispose of the original global list of bodies.
  body_array_free(bodies);

  // Set up a Verlet solver.
  nb->solver = verlet_ode_solver_new(nb->comm, (int)nb->bodies->size,
                                     nb, barnes_hut_accel, NULL);

  // Now create the model.
  model_vtable vtable = {.init = nbody_init,
                         .max_dt = nbody_max_dt,
                         .advance = nbody_advance,
                         .finalize = nbody_finalize,
                         .dtor = nbody_dtor};
  return model_new("Barnes-Hut N-body", nb, vtable, MODEL_MPI);
}

// Probe stuff.

typedef struct
{
  char* body_name;
  nbody_t* model;
  size_t body_index;
} nbody_probe_t;

static void nbody_probe_set_model(void* context, void* model_context)
{
  nbody_probe_t* p = context;
  nbody_t* m = model_context;
  p->model = m;
  if (p->body_name != NULL)
  {
    for (size_t b = 0; b < m->bodies->size; ++b)
    {
      if (strcmp(p->body_name, m->bodies->data[b]->name) == 0)
        p->body_index = b;
    }
  }
}

static void nbody_probe_acquire_x(void* context, real_t t, probe_data_t* data)
{
  nbody_probe_t* p = context;
  if (p->body_index != -1)
  {
    nbody_t* m = p->model;
    data->data[0] = m->U[6*p->body_index];
    data->data[1] = m->U[6*p->body_index+1];
    data->data[2] = m->U[6*p->body_index+2];
  }
}

static void nbody_probe_acquire_v(void* context, real_t t, probe_data_t* data)
{
  nbody_probe_t* p = context;
  if (p->body_index != -1)
  {
    nbody_t* m = p->model;
    data->data[0] = m->U[6*p->body_index+3];
    data->data[1] = m->U[6*p->body_index+4];
    data->data[2] = m->U[6*p->body_index+5];
  }
}

static void nbody_probe_acquire_E(void* context, real_t t, probe_data_t* data)
{
  nbody_probe_t* p = context;
  nbody_t* m = p->model;
  data->data[0] = nbody_total_energy(m);
}

static void nbody_probe_dtor(void* context)
{
  nbody_probe_t* p = context;
  if (p->body_name != NULL)
    string_free(p->body_name);
  polymec_free(p);
}

probe_t* nbody_x_probe_new(const char* body_name)
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "%s position", body_name);
  snprintf(data_name, 128, "x_%s", body_name);
  size_t shape[1] = {3};
  nbody_probe_t* nbody_p = polymec_malloc(sizeof(nbody_probe_t));
  nbody_p->body_name = string_dup(body_name);
  nbody_p->body_index = -1;
  probe_vtable vtable = {.set_model = nbody_probe_set_model,
                         .acquire = nbody_probe_acquire_x,
                         .dtor = nbody_probe_dtor};
  return probe_new(probe_name, data_name, 1, shape, nbody_p, vtable);
}

probe_t* nbody_v_probe_new(const char* body_name)
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "%s velocity", body_name);
  snprintf(data_name, 128, "v_%s", body_name);
  size_t shape[1] = {3};
  nbody_probe_t* nbody_p = polymec_malloc(sizeof(nbody_probe_t));
  nbody_p->body_name = string_dup(body_name);
  nbody_p->body_index = -1;
  probe_vtable vtable = {.set_model = nbody_probe_set_model,
                         .acquire = nbody_probe_acquire_v,
                         .dtor = nbody_probe_dtor};
  return probe_new(probe_name, data_name, 1, shape, nbody_p, vtable);
}

probe_t* nbody_E_probe_new()
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "Total energy");
  snprintf(data_name, 128, "E");
  nbody_probe_t* nbody_p = polymec_malloc(sizeof(nbody_probe_t));
  nbody_p->body_name = NULL;
  nbody_p->body_index = -1;
  probe_vtable vtable = {.set_model = nbody_probe_set_model,
                         .acquire = nbody_probe_acquire_E,
                         .dtor = nbody_probe_dtor};
  return probe_new(probe_name, data_name, 0, NULL, nbody_p, vtable);
}

