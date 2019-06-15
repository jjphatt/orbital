// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "verlet_solver.h"

typedef struct
{
  MPI_Comm comm;

  // Requested time step.
  real_t dt;

  // Context and v-table.
  void* context;
  bool (*compute_dvdt)(void* context, real_t t, nvector_t* u, nvector_t* dvdt);
  void (*dtor)(void* context);

  // Bookkeeping.
  int N;
  nvector_t *a_old, *a_new;
} verlet_t;

static bool verlet_step(void* context, real_t max_dt, real_t* t, nvector_t* u)
{
  START_FUNCTION_TIMER();
  verlet_t* v = context;

  // Compute the acceleration at time t.
  real_t t_old = *t;
  bool status = v->compute_dvdt(v->context, t_old, u, v->a_old);
  if (!status)
  {
    scasm_error("verlet_step: evaluation of acceleration failed at t = %g.\n",
                t_old);
  }

  // Compute x(t + dt) and copy it into U.
  size_t N = nvector_local_size(u);
  real_t ui[N], a_old[N/2];
  nvector_get_local_values(u, ui);
  nvector_get_local_values(v->a_old, a_old);
  int n = (int)N/6;
  for (int i = 0; i < n; ++i)
  {
    real_t x_old = ui[6*i];
    real_t y_old = ui[6*i+1];
    real_t z_old = ui[6*i+2];
    real_t vx_old = ui[6*i+3];
    real_t vy_old = ui[6*i+4];
    real_t vz_old = ui[6*i+5];
    real_t ax_old = a_old[3*i];
    real_t ay_old = a_old[3*i+1];
    real_t az_old = a_old[3*i+2];
    real_t x_new = x_old + vx_old * max_dt + 0.5 * ax_old * max_dt * max_dt;
    real_t y_new = y_old + vy_old * max_dt + 0.5 * ay_old * max_dt * max_dt;
    real_t z_new = z_old + vz_old * max_dt + 0.5 * az_old * max_dt * max_dt;
    ui[6*i] = x_new;
    ui[6*i+1] = y_new;
    ui[6*i+2] = z_new;
  }
  nvector_set_local_values(u, ui);

  // Compute the acceleration at time t + dt.
  real_t t_new = t_old + max_dt;
  status = v->compute_dvdt(v->context, t_new, u, v->a_new);
  if (!status)
  {
    scasm_error("verlet_step: evaluation of acceleration failed at t = %g.\n",
                t_new);
  }

  // Compute v(t + dt) and copy it into U.
  real_t a_new[N/2];
  nvector_get_local_values(v->a_new, a_new);
  for (int i = 0; i < N; ++i)
  {
    real_t vx_old = ui[6*i+3];
    real_t vy_old = ui[6*i+4];
    real_t vz_old = ui[6*i+5];
    real_t ax_old = a_old[3*i];
    real_t ay_old = a_old[3*i+1];
    real_t az_old = a_old[3*i+2];
    real_t ax_new = a_new[3*i];
    real_t ay_new = a_new[3*i+1];
    real_t az_new = a_new[3*i+2];
    real_t vx_new = vx_old + 0.5 * (ax_old + ax_new) * max_dt;
    real_t vy_new = vy_old + 0.5 * (ay_old + ay_new) * max_dt;
    real_t vz_new = vz_old + 0.5 * (az_old + az_new) * max_dt;
    ui[6*i+3] = vx_new;
    ui[6*i+4] = vy_new;
    ui[6*i+5] = vz_new;
  }
  nvector_set_local_values(u, ui);

  // Finish up and get out.
  *t += max_dt;
  STOP_FUNCTION_TIMER();
  return true;
}

static void verlet_dtor(void* context)
{
  verlet_t* v = context;
  scasm_free(v->a_old);
  scasm_free(v->a_new);
  if ((v->context != NULL) && (v->dtor != NULL))
    v->dtor(v->context);
  scasm_free(v);
}

ode_solver_t* verlet_solver_new(MPI_Comm comm, int N, void* context,
                                bool (*compute_dvdt)(void* context, real_t t,
                                                     nvector_t* u, nvector_t* dvdt),
                                void (*dtor)(void* context))
{
  ASSERT(N > 0);
  ASSERT(compute_dvdt != NULL);

  verlet_t* v = scasm_malloc(sizeof(verlet_t));
  v->comm = comm;
  v->N = N;
  v->context = context;
  v->compute_dvdt = compute_dvdt;
  v->dtor = dtor;
  v->a_old = scasm_malloc(sizeof(real_t) * 3 * N);
  v->a_new = scasm_malloc(sizeof(real_t) * 3 * N);

  ode_solver_methods vtable = {.step = verlet_step, .dtor = verlet_dtor};
  int order = 2;
  return ode_solver_new("Verlet ODE solver", v, vtable, order);
}

