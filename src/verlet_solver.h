// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ORBITAL_VERLET_SOLVER_H
#define ORBITAL_VERLET_SOLVER_H

#include "orbital.h"

// This type of ODE solver integrates a Newtonian system
// using the Velocity Verlet method.

// Creates a Verlet solver that integrates a solution vector u = (x, v) by
// computing accelerations dv/dt.
ode_solver_t* verlet_solver_new(void* context, nvector_t* u,
                                void (*compute_dvdt)(void* context, real_t t,
                                                     nvector_t* u, nvector_t* dvdt),
                                void (*dtor)(void* context));

#endif

