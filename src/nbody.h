// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#ifndef ORBITAL_NBODY_H
#define ORBITAL_NBODY_H

#include "body.h"

// This constructs a brute-force N-body gravity model with the given 
// Gravitational constant G and a list of bodies.
model_t* brute_force_nbody_new(real_t G,
                               body_array_t* bodies);

// This constructs a Barnes-Hut N-body gravity model with the given 
// Gravitational constant G, the given far-field parameter theta, and a 
// list of bodies.
model_t* barnes_hut_nbody_new(real_t G,
                              real_t theta,
                              body_array_t* bodies);

// Creates a probe that tracks the position of the body with the given name.
probe_t* nbody_x_probe_new(const char* body_name);

// Creates a probe that tracks the velocity of the body with the given name.
probe_t* nbody_v_probe_new(const char* body_name);

// Creates a probe that tracks the total energy of the system.
probe_t* nbody_E_probe_new(void);

#endif

