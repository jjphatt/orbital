// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "verlet_ode_solver.h"

typedef struct
{
  MPI_Comm comm;

  // Requested time step.
  real_t dt;

  // Context and v-table.
  void* context; 
  int (*compute_dvdt)(void* context, real_t t, real_t* U, real_t* dvdt);
  void (*dtor)(void* context);

  // Bookkeeping.
  int N;
  real_t *a_old, *a_new;
} verlet_t;

static bool verlet_step(void* context, real_t max_dt, real_t* t, real_t* U)
{
  START_FUNCTION_TIMER();
  verlet_t* v = context;

  int N = v->N;

  // Compute the acceleration at time t.
  real_t t_old = *t;
  int status = v->compute_dvdt(v->context, t_old, U, v->a_old);
  if (status != 0)
    polymec_error("verlet_step: evaluation of acceleration failed at t = %g.\n", t_old);

  // Compute x(t + dt) and copy it into U.
  for (int i = 0; i < N; ++i)
  {
    real_t x_old = U[6*i];
    real_t y_old = U[6*i+1];
    real_t z_old = U[6*i+2];
    real_t vx_old = U[6*i+3];
    real_t vy_old = U[6*i+4];
    real_t vz_old = U[6*i+5];
    real_t ax_old = v->a_old[3*i];
    real_t ay_old = v->a_old[3*i+1];
    real_t az_old = v->a_old[3*i+2];
    real_t x_new = x_old + vx_old * max_dt + 0.5 * ax_old * max_dt * max_dt;
    real_t y_new = y_old + vy_old * max_dt + 0.5 * ay_old * max_dt * max_dt;
    real_t z_new = z_old + vz_old * max_dt + 0.5 * az_old * max_dt * max_dt;
    U[6*i] = x_new;
    U[6*i+1] = y_new;
    U[6*i+2] = z_new;
  }

  // Compute the acceleration at time t + dt.
  real_t t_new = t_old + max_dt;
  status = v->compute_dvdt(v->context, t_new, U, v->a_new);
  if (status != 0)
    polymec_error("verlet_step: evaluation of acceleration failed at t = %g.\n", t_new);

  // Compute v(t + dt) and copy it into U.
  for (int i = 0; i < N; ++i)
  {
    real_t vx_old = U[6*i+3];
    real_t vy_old = U[6*i+4];
    real_t vz_old = U[6*i+5];
    real_t ax_old = v->a_old[3*i];
    real_t ay_old = v->a_old[3*i+1];
    real_t az_old = v->a_old[3*i+2];
    real_t ax_new = v->a_new[3*i];
    real_t ay_new = v->a_new[3*i+1];
    real_t az_new = v->a_new[3*i+2];
    real_t vx_new = vx_old + 0.5 * (ax_old + ax_new) * max_dt;
    real_t vy_new = vy_old + 0.5 * (ay_old + ay_new) * max_dt;
    real_t vz_new = vz_old + 0.5 * (az_old + az_new) * max_dt;
    U[6*i+3] = vx_new;
    U[6*i+4] = vy_new;
    U[6*i+5] = vz_new;
  }

  // Finish up and get out.
  *t += max_dt;
  STOP_FUNCTION_TIMER();
  return true;
}

static bool verlet_advance(void* context, real_t t1, real_t t2, real_t* x)
{
  // Well, if you're not going to be picky, we're going to be greedy.
  real_t max_dt = t2 - t1; // Take one honkin' step.
  real_t t = t1;
  return verlet_step(context, max_dt, &t, x);
}

static void verlet_dtor(void* context)
{
  verlet_t* v = context;
  polymec_free(v->a_old);
  polymec_free(v->a_new);
  if ((v->context != NULL) && (v->dtor != NULL))
    v->dtor(v->context);
  polymec_free(v);
}

ode_solver_t* verlet_ode_solver_new(MPI_Comm comm,
                                    int N,
                                    void* context, 
                                    int (*compute_dvdt)(void* context, real_t t, real_t* U, real_t* dUdt),
                                    void (*dtor)(void* context))
{
  ASSERT(N > 0);
  ASSERT(compute_dvdt != NULL);

  verlet_t* v = polymec_malloc(sizeof(verlet_t));
  v->comm = comm;
  v->N = N;
  v->context = context;
  v->compute_dvdt = compute_dvdt;
  v->dtor = dtor;
  v->a_old = polymec_malloc(sizeof(real_t) * 3 * N);
  v->a_new = polymec_malloc(sizeof(real_t) * 3 * N);

  ode_solver_vtable vtable = {.step = verlet_step, 
                              .advance = verlet_advance, 
                              .dtor = verlet_dtor};

  int order = 2;
  return ode_solver_new("Verlet ODE solver", v, vtable, order, 6*N);
}

