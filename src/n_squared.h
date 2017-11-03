// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#ifndef ORBITAL_N_SQUARED_H
#define ORBITAL_N_SQUARED_H

#include "body.h"

// This constructs our "N-squared" gravity model.
model_t* n_squared_new(real_t G,
                       body_array_t* bodies);

// Creates a probe that tracks the position of the body with the given name.
probe_t* n_squared_x_probe_new(const char* body_name);

// Creates a probe that tracks the velocity of the body with the given name.
probe_t* n_squared_v_probe_new(const char* body_name);
#endif

