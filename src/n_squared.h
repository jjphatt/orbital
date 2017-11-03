// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#ifndef ORBITAL_N_SQUARED_H
#define ORBITAL_N_SQUARED_H

#include "body.h"

// This constructs our "N-squared" gravity model.
model_t* n_squared_new(real_t G,
                       body_array_t* bodies);

#endif

