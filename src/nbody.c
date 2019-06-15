// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#include "nbody.h"
#include "verlet_solver.h"
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

  // Time integrator stuff.
  bool (*accel)(void* context, real_t t, nvector_t* U, nvector_t* dvdt);
  ode_solver_t* solver;

  // Solution data.
  nvector_t* u;
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
    real_t vx = b1->v.x;
    real_t vy = b1->v.y;
    real_t vz = b1->v.z;
    E += 0.5 * m1 * (vx*vx + vy*vy + vz*vz);
    point_t x1 = b1->x;
    for (size_t bb = 0; bb < nb->bodies->size; ++bb)
    {
      if (b == bb) continue;
      body_t* b2 = nb->bodies->data[bb];
      real_t m2 = b2->m;
      point_t x2 = b2->x;
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

  // Count up all the bodies.
  int local_size = (int)(6*nb->bodies->size);
  index_t global_size = local_size;
  MPI_Allreduce(MPI_IN_PLACE, &global_size, 1, MPI_INDEX_T, MPI_SUM, nb->comm);

  // Create and populate the solution vector.
  if (nb->u == NULL)
    nb->u = nvector_new(nb->comm, local_size, global_size);
  real_t u[local_size];
  nb->v_min = REAL_MAX;
  for (size_t b = 0; b < nb->bodies->size; ++b)
  {
    body_t* body = nb->bodies->data[b];
    u[6*b]   = body->x.x;
    u[6*b+1] = body->x.y;
    u[6*b+2] = body->x.z;
    u[6*b+3] = body->v.x;
    u[6*b+4] = body->v.y;
    u[6*b+5] = body->v.z;
    nb->v_min = MIN(nb->v_min, vector_mag(&body->v));
  }
  nvector_set_local_values(nb->u, u);

  // Compute the initial energy of the system.
  nb->E0 = nbody_total_energy(nb);

  // Set up a Verlet solver.
  if (nb->solver != NULL)
    ode_solver_free(nb->solver);
  nb->solver = verlet_solver_new(nb->comm, (int)nb->bodies->size,
                                 nb, nb->accel, NULL);
}

static real_t nbody_max_dt(void* context, real_t t, char* reason)
{
//  nbody_t* nb = context;

  // The timestep is limited by the ratio of the maximum acceleration to
  // the minimum speed.
  sprintf(reason, "ratio of max acceleration to minimum speed");
  return REAL_MAX;//0.2 * (nb->a_max / nb->v_min);
}

static real_t nbody_advance(void* context, real_t max_dt, real_t t)
{
  nbody_t* nb = context;
  scasm_suspend_fpe();
  bool solved = ode_solver_advance(nb->solver, t, t + max_dt, nb->u);
  scasm_restore_fpe();
  if (!solved)
    return 0.0;
  else
  {
    // Copy the data from the solution vector into the bodies.
    real_t u[6*nb->bodies->size];
    nvector_get_local_values(nb->u, u);
    for (size_t b = 0; b < nb->bodies->size; ++b)
    {
      body_t* body = nb->bodies->data[b];
      body->x.x = u[6*b];
      body->x.y = u[6*b+1];
      body->x.z = u[6*b+2];
      body->v.x = u[6*b+3];
      body->v.y = u[6*b+4];
      body->v.z = u[6*b+5];
      nb->v_min = MIN(nb->v_min, vector_mag(&body->v));
    }
  }
  return max_dt;
}

static void nbody_finalize(void* context, int step, real_t t)
{
  // Compute the final energy and compare it to the initial energy.
  nbody_t* nb = context;
  real_t E = nbody_total_energy(nb);
  log_detail("Fractional energy change: %g", (E - nb->E0) / nb->E0 - 1.0);
}

static void nbody_plot(void* context,
                       const char* file_prefix,
                       const char* directory,
                       real_t time,
                       int step)
{
  START_FUNCTION_TIMER();

  nbody_t* nb = context;

  // Open a SILO viz file.
  silo_file_t* silo = silo_file_new(nb->comm, file_prefix, directory,
                                    1, step, time);

  // Write a point cloud to the file.
  int N = (int)(nb->bodies->size);
  point_cloud_t* cloud = point_cloud_new(nb->comm, N);
  for (int i = 0; i < N; ++i)
    cloud->points[i] = nb->bodies->data[i]->x;
  silo_file_write_point_cloud(silo, "bodies", cloud);

  // Write the masses and velocities.
  point_cloud_field_t* field = point_cloud_field_new(cloud, 4);
  DECLARE_POINT_CLOUD_FIELD_ARRAY(f, field);
  for (int i = 0; i < N; ++i)
  {
    body_t* b = nb->bodies->data[i];
    f[i][0] = b->m;
    f[i][1] = b->v.x;
    f[i][2] = b->v.y;
    f[i][3] = b->v.z;
  }
  field_metadata_t* md = point_cloud_field_metadata(field);
  field_metadata_set_name(md, 0, "masses");
  field_metadata_set_name(md, 1, "vx");
  field_metadata_set_name(md, 2, "vy");
  field_metadata_set_name(md, 3, "vz");
  field_metadata_set_vector(md, 1);
  silo_file_write_point_field(silo, "data", "bodies", field);
  silo_file_write_vector_expression(silo, "velocities", "{vx, vy, vz}");

  // Write the file.
  silo_file_close(silo);

  // Clean up.
  point_cloud_field_free(field);
  point_cloud_free(cloud);

  STOP_FUNCTION_TIMER();
}

static void nbody_save(void* context,
                       const char* file_prefix,
                       const char* directory,
                       real_t time,
                       int step)
{
  START_FUNCTION_TIMER();
  nbody_t* nb = context;

  // Open a SILO save file.
  silo_file_t* silo = silo_file_new(nb->comm, file_prefix, directory,
                                    1, step, time);

  // Write the solution vector directly to the file.
  real_t u[6*nb->bodies->size];
  nvector_get_local_values(nb->u, u);
  silo_file_write_real_array(silo, "u", u, 6*nb->bodies->size);

  // Write the masses and body names.
  int N = (int)(nb->bodies->size);
  real_t* masses = scasm_malloc(sizeof(real_t) * N);
  char_array_t* all_names = char_array_new();
  for (int b = 0; b < N; ++b)
  {
    masses[b] = nb->bodies->data[b]->m;
    size_t name_len = strlen(nb->bodies->data[b]->name);
    for (size_t c = 0; c < name_len; ++c)
      char_array_append(all_names, nb->bodies->data[b]->name[c]);
    char_array_append(all_names, ';');
  }
  char_array_append(all_names, '\0');
  silo_file_write_real_array(silo, "m", masses, 6*nb->bodies->size);
  silo_file_write_string(silo, "names", all_names->data);
  silo_file_write_real_array(silo, "E0", &(nb->E0), 1);

  // Write the file.
  silo_file_close(silo);

  // Clean up.
  scasm_free(masses);
  char_array_free(all_names);

  STOP_FUNCTION_TIMER();
}

static bool nbody_load(void* context, const char* file_prefix,
                       const char* directory, real_t* time, int step)
{
  START_FUNCTION_TIMER();
  nbody_t* nb = context;

  // Open a SILO save file.
  silo_file_t* silo = silo_file_open(nb->comm, file_prefix, directory,
                                     step, time);
  if (silo == NULL)
    return false;

  // Dispose of the bodies.
  if (nb->bodies != NULL)
    body_array_free(nb->bodies);
  nb->bodies = body_array_new();

  // Read the solution vector directly from the file.
  size_t six_N;
  real_t* u = silo_file_read_real_array(silo, "u", &six_N);
  int N = (int)(six_N / 6);

  // Read the masses and body names.
  real_t* masses = silo_file_read_real_array(silo, "m", &six_N);
  ASSERT(six_N == (size_t)(6*N));
  const char* all_names = silo_file_read_string(silo, "names");
  const char* name_ptr = all_names;
  for (int b = 0; b < N; ++b)
  {
    // Retrieve the body's name.
    char* semicolon = strstr(name_ptr, ";");
    size_t name_len = semicolon - name_ptr;
    char name[name_len+1];
    memcpy(name, name_ptr, sizeof(char) * name_len);
    name[name_len] = '\0';
    name_ptr = semicolon + 1;

    // Get the other parameters.
    real_t m = masses[b];
    point_t x = {.x = u[6*b], .y = u[6*b+1], .z = u[6*b+2]};
    vector_t v = {.x = u[6*b+3], .y = u[6*b+4], .z = u[6*b+5]};
    body_array_append(nb->bodies, body_new(name, m, &x, &v));
  }

  // Retrieve the initial energy of the system.
  size_t one;
  real_t* E0 = silo_file_read_real_array(silo, "E0", &one);
  ASSERT(one == 1);
  nb->E0 = *E0;

  // Close the file and clean up.
  silo_file_close(silo);

  // Count up all the bodies.
  int local_size = (int)(6*nb->bodies->size);
  index_t global_size = local_size;
  MPI_Allreduce(MPI_IN_PLACE, &global_size, 1, MPI_INDEX_T, MPI_SUM, nb->comm);

  // Create and populate the solution vector.
  if (nb->u == NULL)
    nb->u = nvector_new(nb->comm, local_size, global_size);
  nvector_set_local_values(nb->u, u);
  scasm_free(u);

  // Clean up.
  scasm_free(E0);
  scasm_free(masses);
  string_free(all_names);

  // Set up a Verlet solver.
  if (nb->solver != NULL)
    ode_solver_free(nb->solver);
  nb->solver = verlet_solver_new(nb->comm, (int)nb->bodies->size, nb, nb->accel,
                                 NULL);

  STOP_FUNCTION_TIMER();
  return true;
}

static void nbody_dtor(void* context)
{
  nbody_t* nb = context;
  if (nb->bodies != NULL)
    body_array_free(nb->bodies);
  if (nb->u != NULL)
    nvector_free(nb->u);
  if (nb->solver != NULL)
    ode_solver_free(nb->solver);
  if (nb->tree != NULL)
    barnes_hut_tree_free(nb->tree);
  scasm_free(nb);
}

static body_array_t* partition_bodies(nbody_t* nb, body_array_t* bodies)
{
  START_FUNCTION_TIMER();

  int N = (int)(bodies->size);
  if (N == 0)
  {
    STOP_FUNCTION_TIMER();
    return body_array_new();
  }

  // Partition the bodies and their data. Since each process has all of the
  // data, we just compute a partition vector.
  point_cloud_t* cloud = point_cloud_new(nb->comm, N);
  for (int b = 0; b < N; ++b)
    cloud->points[b] = bodies->data[b]->x;
  int64_t* P = partition_points(cloud->points, cloud->num_points,
                                nb->comm, NULL, 1.05, true);

  // Use the partition vector to cull the off-process bodies.
  body_array_t* local_bodies = body_array_new();
  for (int i = 0; i < N; ++i)
  {
    if (P[i] == nb->rank)
      body_array_append(local_bodies, body_clone(bodies->data[i]));
  }
  nb->bodies = local_bodies;

  // Clean up.
  scasm_free(P);
  point_cloud_free(cloud);

  STOP_FUNCTION_TIMER();
  return local_bodies;
}

// Solver functions.
static bool brute_force_accel(void* context, real_t t, nvector_t* u, nvector_t* dvdt)
{
  nbody_t* nb = context;
  real_t G = nb->G;
  nb->a_max = 0.0;
  nb->v_min = REAL_MAX;

  // Extract point coordinates and masses from solution data.
  int N = (int)(nb->bodies->size);
  point_t x[N];
  real_t m[N], ui[6*N];
  nvector_get_local_values(u, ui);
  for (int i = 0; i < N; ++i)
  {
    x[i].x = ui[6*i];
    x[i].y = ui[6*i+1];
    x[i].z = ui[6*i+2];
    m[i] = nb->bodies->data[i]->m;
  }

  // How many points on other processes?
  int Np[nb->nprocs];
  MPI_Allgather(&N, 1, MPI_INT, Np, 1, MPI_INT, nb->comm);
  int Ntot = 0;
  for (int p = 0; p < nb->nprocs; ++p)
    Ntot += Np[p];

  // Get point data from other processes.
  point_t* all_points = scasm_malloc(sizeof(point_t) * Ntot);
  MPI_Allgather(x, 3*N, MPI_REAL_T, all_points, 3*N, MPI_REAL_T, nb->comm);

  // Now compute forces on all pairs, with i ranging over the entirety of
  // the points, and j over just the local ones.
  vector_t* all_accels = scasm_malloc(sizeof(vector_t) * Ntot);
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
  MPI_Allreduce(MPI_IN_PLACE, all_accels, 3*Ntot, MPI_REAL_T, MPI_SUM, nb->comm);

  // Extract acceleration data for this process.
  int start = 0;
  for (int p = 0; p < nb->rank; ++p)
    start += Np[p];
  real_t a[3*N];
  memcpy(a, &(all_accels[start]), sizeof(vector_t) * N);
  nvector_set_local_values(dvdt, a);

  // Clean up.
  scasm_free(all_accels);
  scasm_free(all_points);

  return 0;
}

static bool barnes_hut_accel(void* context, real_t t, nvector_t* u, nvector_t* dvdt)
{
  nbody_t* nb = context;
  int N = (int)(nb->bodies->size);
  point_t x[N];
  real_t m[N], ui[6*N];
  nvector_get_local_values(u, ui);
  for (int i = 0; i < N; ++i)
  {
    x[i].x = ui[6*i];
    x[i].y = ui[6*i+1];
    x[i].z = ui[6*i+2];
    m[i] = nb->bodies->data[i]->m;
  }
  real_t a[3*N];
  barnes_hut_tree_compute_forces(nb->tree, nb->G, x, m, N, (vector_t*)a);
  for (int i = 0; i < N; ++i)
  {
    a[3*i] /= m[i];
    a[3*i+1] /= m[i];
    a[3*i+2] /= m[i];
  }
  nvector_set_local_values(dvdt, a);

  return 0;
}

//------------------------------------------------------------------------

model_t* brute_force_nbody_new(real_t G,
                               body_array_t* bodies)
{
  ASSERT(G >= 0.0);

  nbody_t* nb = scasm_malloc(sizeof(nbody_t));
  nb->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(nb->comm, &(nb->rank));
  MPI_Comm_size(nb->comm, &(nb->nprocs));
  nb->G = G;
  nb->bodies = partition_bodies(nb, bodies);
  nb->tree = NULL;
  nb->u = NULL;
  nb->solver = NULL;
  nb->accel = brute_force_accel;

  // Dispose of the original global list of bodies.
  body_array_free(bodies);

  // Now create the model.
  model_methods vtable = {.init = nbody_init,
                         .max_dt = nbody_max_dt,
                         .advance = nbody_advance,
                         .finalize = nbody_finalize,
                         .dtor = nbody_dtor};
  return model_new("Brute-Force N-body", nb, vtable, MODEL_MPI);
}

model_t* barnes_hut_nbody_new(real_t G,
                              real_t theta,
                              body_array_t* bodies)
{
  ASSERT(G >= 0.0);
  ASSERT(theta >= 0.0);

  nbody_t* nb = scasm_malloc(sizeof(nbody_t));
  nb->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(nb->comm, &(nb->rank));
  MPI_Comm_size(nb->comm, &(nb->nprocs));
  nb->G = G;
  nb->bodies = partition_bodies(nb, bodies);
  nb->tree = barnes_hut_tree_new(nb->comm, theta);
  nb->u = NULL;
  nb->solver = NULL;
  nb->accel = barnes_hut_accel;

  // Dispose of the original global list of bodies.
  body_array_free(bodies);

  // Now create the model.
  model_methods vtable = {.init = nbody_init,
                         .max_dt = nbody_max_dt,
                         .advance = nbody_advance,
                         .finalize = nbody_finalize,
                         .plot = nbody_plot,
                         .save = nbody_save,
                         .load = nbody_load,
                         .dtor = nbody_dtor};
  return model_new("Barnes-Hut N-body", nb, vtable, MODEL_MPI);
}

// Probe stuff.

typedef struct
{
  const char* body_name;
  nbody_t* model;
  int body_index, body_rank;
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
        p->body_index = (int)b;
    }
  }

  // Find the rank on which the body lives.
  int val_and_rank[2] = {p->body_index, m->rank};
  MPI_Allreduce(MPI_IN_PLACE, val_and_rank, 1, MPI_2INT, MPI_MAXLOC, m->comm);
  p->body_rank = val_and_rank[1];
}

static void nbody_probe_acquire_x(void* context, real_t t, probe_data_t* data)
{
  nbody_probe_t* p = context;
  nbody_t* m = p->model;
  if (p->body_index != -1)
  {
    body_t* b = m->bodies->data[p->body_index];
    data->data[0] = b->x.x;
    data->data[1] = b->x.y;
    data->data[2] = b->x.z;
  }
  MPI_Bcast(data->data, 3, MPI_REAL_T, p->body_rank, m->comm);
}

static void nbody_probe_acquire_v(void* context, real_t t, probe_data_t* data)
{
  nbody_probe_t* p = context;
  nbody_t* m = p->model;
  if (p->body_index != -1)
  {
    body_t* b = m->bodies->data[p->body_index];
    data->data[0] = b->v.x;
    data->data[1] = b->v.y;
    data->data[2] = b->v.z;
  }
  MPI_Bcast(data->data, 3, MPI_REAL_T, p->body_rank, m->comm);
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
  scasm_free(p);
}

probe_t* nbody_x_probe_new(const char* body_name)
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "%s position", body_name);
  snprintf(data_name, 128, "x_%s", body_name);
  size_t shape[1] = {3};
  nbody_probe_t* nbody_p = scasm_malloc(sizeof(nbody_probe_t));
  nbody_p->body_name = string_dup(body_name);
  nbody_p->body_index = -1;
  probe_methods vtable = {.set_model = nbody_probe_set_model,
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
  nbody_probe_t* nbody_p = scasm_malloc(sizeof(nbody_probe_t));
  nbody_p->body_name = string_dup(body_name);
  nbody_p->body_index = -1;
  probe_methods vtable = {.set_model = nbody_probe_set_model,
                         .acquire = nbody_probe_acquire_v,
                         .dtor = nbody_probe_dtor};
  return probe_new(probe_name, data_name, 1, shape, nbody_p, vtable);
}

probe_t* nbody_E_probe_new()
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "Total energy");
  snprintf(data_name, 128, "E");
  nbody_probe_t* nbody_p = scasm_malloc(sizeof(nbody_probe_t));
  nbody_p->body_name = NULL;
  nbody_p->body_index = -1;
  probe_methods vtable = {.set_model = nbody_probe_set_model,
                          .acquire = nbody_probe_acquire_E,
                          .dtor = nbody_probe_dtor};
  return probe_new(probe_name, data_name, 0, NULL, nbody_p, vtable);
}

