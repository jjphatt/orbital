// Copyright (c) 2016, Jeffrey N. Johnson
// All rights reserved.

#ifndef ORBITAL_BODY_H
#define ORBITAL_BODY_H

#include "polymec.h"

// A body floats through 3D space under the influence of gravity.
// Objects of this type are garbage-collected.
typedef struct 
{
  char* name;

  real_t m;   // mass
  point_t x;  // position
  vector_t v; // velocity
} body_t;

// Creates a new body with a descriptive name, a mass m, a position x, 
// and a velocity v.
body_t* body_new(const char* name,
                 real_t m,
                 point_t* x, 
                 vector_t* v);

DEFINE_ARRAY(body_array, body_t*)

#endif

