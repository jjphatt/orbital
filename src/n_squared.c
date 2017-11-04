// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#include "n_squared.h"

typedef struct
{
  real_t G;
  body_array_t* bodies;
  ode_solver_t* solver;
  real_t* U;
  real_t v_min, a_max;
} nsq_t;

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

  // Compute the maximum acceleration on any body.
  real_t G = nsq->G;
  nsq->a_max = 0.0;
  for (size_t b = 0; b < nsq->bodies->size; ++b)
  {
//    body_t* b1 = nsq->bodies->data[b];
    point_t x1 = {.x = nsq->U[6*b], .y = nsq->U[6*b+1], .z = nsq->U[6*b+2]};
    for (size_t bb = 0; bb < nsq->bodies->size; ++bb)
    {
      if (b == bb) continue;
      body_t* b2 = nsq->bodies->data[bb];
      point_t x2 = {.x = nsq->U[6*bb], .y = nsq->U[6*bb+1], .z = nsq->U[6*bb+2]};
      vector_t r12;
      point_displacement(&x1, &x2, &r12);
      real_t r2 = vector_dot(&r12, &r12);
      nsq->a_max = MAX(nsq->a_max, G * b2->m / r2);
    }
  }

  if (nsq->v_min == 0.0)
    nsq->v_min = 0.2 * nsq->a_max;
}

static real_t nsq_max_dt(void* context, real_t t, char* reason)
{
  nsq_t* nsq = context;

  // The timestep is limited by the ratio of the maximum acceleration to 
  // the minimum speed.
  return 0.2 * (nsq->a_max / nsq->v_min);
}

static real_t nsq_advance(void* context, real_t max_dt, real_t t)
{
//  nsq_t* nsq = context;
  return max_dt;
}

static void nsq_finalize(void* context, int step, real_t t)
{
}

static void nsq_add_probe(void* context, void* probe_context)
{
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
static int nsq_accel(void* context, real_t t, real_t* U, real_t* dUdt)
{
  nsq_t* nsq = context;

  real_t G = nsq->G;
  nsq->a_max = 0.0;
  nsq->v_min = REAL_MAX;
  int N = (int)(nsq->bodies->size);
  for (size_t i = 0; i < N; ++i)
  {
    point_t x1 = {.x = U[6*i], .y = U[6*i+1], .z = U[6*i+2]};
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
      a.x -= G * m2 * r12.x * r3_inv;
      a.y -= G * m2 * r12.y * r3_inv;
      a.z -= G * m2 * r12.z * r3_inv;
      nsq->a_max = MAX(nsq->a_max, vector_mag(&a));
    }

    // Compute derivatives.
    dUdt[6*i]   = U[6*i+3];
    dUdt[6*i+1] = U[6*i+4];
    dUdt[6*i+2] = U[6*i+5];
    dUdt[6*i+3] = a.x;
    dUdt[6*i+4] = a.y;
    dUdt[6*i+5] = a.z;
  }
  return 0;
}

static real_t nsq_stable_dt(void* context, real_t t, real_t* U)
{
  nsq_t* nsq = context;
  return 0.2 * (nsq->a_max / nsq->v_min);
}

// API
model_t* n_squared_new(real_t G,
                       body_array_t* bodies)
{
  ASSERT(G >= 0.0);

  nsq_t* nsq = polymec_malloc(sizeof(nsq_t));
  nsq->G = G;
  nsq->bodies = bodies;
  nsq->U = polymec_malloc(sizeof(real_t) * 6 * bodies->size);

  // Set up an explicit RK4 solver.
  nsq->solver = explicit_ark_ode_solver_new(4, MPI_COMM_SELF, 
                                            (int)(6 * bodies->size), 0,
                                            nsq, nsq_accel, nsq_stable_dt, 
                                            NULL);

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

  // Now create the model.
  model_vtable vtable = {.init = nsq_init,
                         .max_dt = nsq_max_dt,
                         .advance = nsq_advance,
                         .finalize = nsq_finalize,
                         .add_probe = nsq_add_probe,
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
  for (size_t b = 0; b < m->bodies->size; ++b)
  {
    if (strcmp(p->body_name, m->bodies->data[b]->name) == 0)
      p->body_index = b;
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

static void nsq_probe_dtor(void* context)
{
  nsq_probe_t* p = context;
  string_free(p->body_name);
  polymec_free(p);
}

probe_t* n_squared_x_probe_new(const char* body_name)
{
  char probe_name[129], data_name[129];
  snprintf(probe_name, 128, "%s position", body_name);
  snprintf(data_name, 128, "%s_x", body_name);
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
  snprintf(data_name, 128, "%s_v", body_name);
  size_t shape[1] = {3};
  nsq_probe_t* nsq_p = polymec_malloc(sizeof(nsq_probe_t));
  nsq_p->body_name = string_dup(body_name);
  nsq_p->body_index = -1;
  probe_vtable vtable = {.set_model = nsq_probe_set_model,
                         .acquire = nsq_probe_acquire_v,
                         .dtor = nsq_probe_dtor};
  return probe_new(probe_name, data_name, 1, shape, nsq_p, vtable);
}

