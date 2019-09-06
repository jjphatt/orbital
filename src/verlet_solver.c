// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "verlet_solver.h"

typedef struct
{
  // Context and v-table.
  void* context;
  void (*compute_dvdt)(void* context, real_t t, nvector_t* u, nvector_t* dvdt);
  void (*dtor)(void* context);

  // Bookkeeping.
  nvector_t *a_old, *a_new;
} verlet_t;

static bool verlet_step(void* context, real_t max_dt, real_t* t, nvector_t* u)
{
  START_FUNCTION_TIMER();
  verlet_t* v = context;

  // Compute the acceleration at time t.
  real_t t_old = *t;
  v->compute_dvdt(v->context, t_old, u, v->a_old);

  // Compute x(t + dt) and copy it into u.
  int N = parallel_real_nvector_local_size(u);
  real_t* ui = parallel_real_nvector_local_data(u);
  real_t* a_old = parallel_real_nvector_local_data(v->a_old);
  real_t* a_new = parallel_real_nvector_local_data(v->a_new);
  int n = N/6;
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

  // Compute the acceleration at time t + dt.
  real_t t_new = t_old + max_dt;
  v->compute_dvdt(v->context, t_new, u, v->a_new);

  // Compute v(t + dt) and copy it into u.
  for (int i = 0; i < n; ++i)
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

  // Finish up and get out.
  *t += max_dt;
  STOP_FUNCTION_TIMER();
  return true;
}

static void verlet_dtor(void* context)
{
  verlet_t* v = context;
  nvector_free(v->a_old);
  nvector_free(v->a_new);
  if ((v->context != NULL) && (v->dtor != NULL))
    v->dtor(v->context);
  scasm_free(v);
}

ode_solver_t* verlet_solver_new(void* context, nvector_t* u,
                                void (*compute_dvdt)(void* context, real_t t,
                                                     nvector_t* u, nvector_t* dvdt),
                                void (*dtor)(void* context))
{
  ASSERT(u != NULL);
  ASSERT(compute_dvdt != NULL);

  verlet_t* v = scasm_malloc(sizeof(verlet_t));
  v->context = context;
  v->compute_dvdt = compute_dvdt;
  v->dtor = dtor;

  // Set up acceleration n-vectors. These are half the size of u, since they
  // only store accelerations (dv/dt, and not dx/dt).
  index_t global_size = nvector_dimension(u);
  MPI_Comm comm = parallel_real_nvector_comm(u);
  int local_size = parallel_real_nvector_local_size(u);
  v->a_old = parallel_real_nvector_new(comm, local_size/2, global_size/2);
  v->a_new = nvector_clone(v->a_old);

  ode_solver_methods vtable = {.step = verlet_step, .dtor = verlet_dtor};
  int order = 2;
  return ode_solver_new("Verlet ODE solver", v, vtable, order);
}

