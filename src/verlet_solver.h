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

// Creates a Verlet solver that integrates a solution vector representing a
// set of N points, each defined by a position and a velocity in 3D space.
// The ith point is defined by
//   * x, y, and z positions stored in U[6*i], U[6*i+1], U[6*i+2]
//   * x, y, and z velocity components stored in U[6*i+3], U[6*i+4], U[6*i+5]
// Below, compute_dvdt computes the 3 components of the acceleration, placing
// them into dvdt.
ode_solver_t* verlet_solver_new(MPI_Comm comm, int N, void* context,
                                bool (*compute_dvdt)(void* context, real_t t,
                                                     nvector_t* u, nvector_t* dvdt),
                                void (*dtor)(void* context));

#endif

