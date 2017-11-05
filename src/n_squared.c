// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#include "n_squared.h"
#include "verlet_ode_solver.h"

typedef struct
{
  real_t G;
  body_array_t* bodies;
  ode_solver_t* solver;
  real_t* U;
  real_t v_min, a_max;
  real_t E0; // initial energy.
} nsq_t;

static real_t nsq_total_energy(nsq_t* nsq)
{
  // The total energy in the system is the sum of gravitational 
  // potentials plus kinetic energies.
  real_t E = 0.0, G = nsq->G;
  for (size_t b = 0; b < nsq->bodies->size; ++b)
  {
    body_t* b1 = nsq->bodies->data[b];
    real_t m1 = b1->m;
    real_t vx = nsq->U[6*b+3];
    real_t vy = nsq->U[6*b+4];
    real_t vz = nsq->U[6*b+5];
    E += 0.5 * m1 * (vx*vx + vy*vy + vz*vz);
    point_t x1 = {.x = nsq->U[6*b], .y = nsq->U[6*b+1], .z = nsq->U[6*b+2]};
    for (size_t bb = 0; bb < nsq->bodies->size; ++bb)
    {
      if (b == bb) continue;
      body_t* b2 = nsq->bodies->data[bb];
      real_t m2 = b2->m;
      point_t x2 = {.x = nsq->U[6*bb], .y = nsq->U[6*bb+1], .z = nsq->U[6*bb+2]};
      real_t r2 = point_square_distance(&x1, &x2);
      E += G * m1 * m2 / r2;
    }
  }
  return E;
}

static void nsq_init(void* context, real_t t)
{
  nsq_t* nsq = context;
  nsq->v_min = REAL_MAX;
  for (size_t b = 0; b < nsq->bodies->size; ++b)
  {
    body_t* body = nsq->bodies->data[b];
    nsq->U[6*b]   = body->x.x;
    nsq->U[6*b+1] = body->x.y;
    nsq->U[6*b+2] = body->x.z;
    nsq->U[6*b+3] = body->v.x;
    nsq->U[6*b+4] = body->v.y;
    nsq->U[6*b+5] = body->v.z;
    nsq->v_min = MIN(nsq->v_min, vector_mag(&body->v));
  }

  // Compute the initial energy of the system.
  nsq->E0 = nsq_total_energy(nsq);
}

static real_t nsq_max_dt(void* context, real_t t, char* reason)
{
//  nsq_t* nsq = context;

  // The timestep is limited by the ratio of the maximum acceleration to 
  // the minimum speed.
  snprintf(reason, POLYMEC_MODEL_MAXDT_REASON_SIZE,
           "ratio of max acceleration to minimum speed");
  return REAL_MAX;//0.2 * (nsq->a_max / nsq->v_min);
}

static real_t nsq_advance(void* context, real_t max_dt, real_t t)
{
  nsq_t* nsq = context;
  polymec_suspend_fpe();
  bool solved = ode_solver_advance(nsq->solver, t, t + max_dt, nsq->U);
  polymec_restore_fpe();
  if (!solved)
    return 0.0;
  else
    return max_dt;
}

static void nsq_finalize(void* context, int step, real_t t)
{
  // Compute the final energy and compare it to the initial energy.
  nsq_t* nsq = context;
  real_t E = nsq_total_energy(nsq);
  log_detail("Fractional energy change: %g", (E - nsq->E0) / nsq->E0 - 1.0);
}

static void nsq_dtor(void* context)
{
  nsq_t* nsq = context;
  body_array_free(nsq->bodies);
  polymec_free(nsq->U);
  ode_solver_free(nsq->solver);
  polymec_free(nsq);
}

// Solver functions.
static int nsq_accel(void* context, real_t t, real_t* U, real_t* dvdt)
{
  nsq_t* nsq = context;
  real_t G = nsq->G;
  nsq->a_max = 0.0;
  nsq->v_min = REAL_MAX;
  int N = (int)(nsq->bodies->size);
  for (size_t i = 0; i < N; ++i)
  {
    point_t x1 = {.x = U[6*i], .y = U[6*i+1], .z = U[6*i+2]};
    vector_t v1 = {.x = U[6*i+3], .y = U[6*i+4], .z = U[6*i+5]};
    vector_t a = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (size_t j = 0; j < N; ++j)
    {
      if (i == j) continue;
      real_t m2 = nsq->bodies->data[j]->m;
      point_t x2 = {.x = U[6*j], .y = U[6*j+1], .z = U[6*j+2]};
      vector_t r12;
      point_displacement(&x1, &x2, &r12);
      real_t r = vector_mag(&r12);
      real_t r3_inv = 1.0 / (r*r*r);
      a.x += G * m2 * r12.x * r3_inv;
      a.y += G * m2 * r12.y * r3_inv;
      a.z += G * m2 * r12.z * r3_inv;
    }
    nsq->a_max = MAX(nsq->a_max, vector_mag(&a));
    nsq->v_min = MIN(nsq->v_min, vector_mag(&v1));

    // Update the acceleration.
    dvdt[3*i]   = a.x;
    dvdt[3*i+1] = a.y;
    dvdt[3*i+2] = a.z;
  }
//  log_debug("a_max = %g\n", nsq->a_max);
  return 0;
}

//------------------------------------------------------------------------

model_t* n_squared_new(real_t G,
                       body_array_t* bodies)
{
  ASSERT(G >= 0.0);

  // Find out whether we have a Schwartzchild body.
  int sc_index = -1;
  for (size_t b = 0; b < bodies->size; ++b)
  {
    body_t* body = bodies->data[b];
    if (body->schwartzchild)
    {
      ASSERT(sc_index == -1); // Only one of these allowed!
      sc_index = (int)b;
    }
  }

  // Swap any Schwartzchild body with the first one, to impose an ordering.
  if (sc_index != -1)
    body_array_swap(bodies, 0, sc_index);

  nsq_t* nsq = polymec_malloc(sizeof(nsq_t));
  nsq->G = G;
  nsq->bodies = bodies;
  nsq->U = polymec_malloc(sizeof(real_t) * 6 * bodies->size);

  // Set up a Verlet solver.
  nsq->solver = verlet_ode_solver_new(MPI_COMM_SELF, (int)bodies->size,
                                      nsq, nsq_accel, NULL);

  // Now create the model.
  model_vtable vtable = {.init = nsq_init,
                         .max_dt = nsq_max_dt,
                         .advance = nsq_advance,
                         .finalize = nsq_finalize,
                         .dtor = nsq_dtor};
  return model_new("N-squared", nsq, vtable, MODEL_SERIAL);
}

// Probe stuff.

typedef struct
{
  char* body_name;
  nsq_t* model;
  size_t body_index;
} nsq_probe_t;

static void nsq_probe_set_model(void* context, void* model_context)
{
  nsq_probe_t* p = context;
  nsq_t* m = model_context;
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

static void nsq_probe_acquire_x(void* context, real_t t, probe_data_t* data)
{
  nsq_probe_t* p = context;
  if (p->body_index != -1)
  {
    nsq_t* m = p->model;
    data->data[0] = m->U[6*p->body_index];
    data->data[1] = m->U[6*p->body_index+1];
    data->data[2] = m->U[6*p->body_index+2];
  }
}

static void nsq_probe_acquire_v(void* context, real_t t, probe_data_t* data)
{
  nsq_probe_t* p = context;
  if (p->body_index != -1)
  {
    nsq_t* m = p->model;
    data->data[0] = m->U[6*p->body_index+3];
    data->data[1] = m->U[6*p->body_index+4];
    data->data[2] = m->U[6*p->body_index+5];
  }
}

static void nsq_probe_acquire_E(void* context, real_t t, probe_data_t* data)
{
  nsq_probe_t* p = context;
  nsq_t* m = p->model;
  data->data[0] = nsq_total_energy(m);
}

static void nsq_probe_dtor(void* context)
{
  nsq_probe_t* p = context;
  if (p->body_name != NULL)
    string_free(p->body_name);
  polymec_free(p);
}

probe_t* n_squared_x_probe_new(const char* body_name)
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "%s position", body_name);
  snprintf(data_name, 128, "x_%s", body_name);
  size_t shape[1] = {3};
  nsq_probe_t* nsq_p = polymec_malloc(sizeof(nsq_probe_t));
  nsq_p->body_name = string_dup(body_name);
  nsq_p->body_index = -1;
  probe_vtable vtable = {.set_model = nsq_probe_set_model,
                         .acquire = nsq_probe_acquire_x,
                         .dtor = nsq_probe_dtor};
  return probe_new(probe_name, data_name, 1, shape, nsq_p, vtable);
}

probe_t* n_squared_v_probe_new(const char* body_name)
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "%s velocity", body_name);
  snprintf(data_name, 128, "v_%s", body_name);
  size_t shape[1] = {3};
  nsq_probe_t* nsq_p = polymec_malloc(sizeof(nsq_probe_t));
  nsq_p->body_name = string_dup(body_name);
  nsq_p->body_index = -1;
  probe_vtable vtable = {.set_model = nsq_probe_set_model,
                         .acquire = nsq_probe_acquire_v,
                         .dtor = nsq_probe_dtor};
  return probe_new(probe_name, data_name, 1, shape, nsq_p, vtable);
}

probe_t* n_squared_E_probe_new()
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "Total energy");
  snprintf(data_name, 128, "E");
  nsq_probe_t* nsq_p = polymec_malloc(sizeof(nsq_probe_t));
  nsq_p->body_name = NULL;
  nsq_p->body_index = -1;
  probe_vtable vtable = {.set_model = nsq_probe_set_model,
                         .acquire = nsq_probe_acquire_E,
                         .dtor = nsq_probe_dtor};
  return probe_new(probe_name, data_name, 0, NULL, nsq_p, vtable);
}

