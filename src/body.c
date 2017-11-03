// Copyright (c) 2016, Jeffrey N. Johnson
// All rights reserved.

#include "body.h"

body_t* body_new(const char* name,
                 real_t m,
                 point_t* x, 
                 vector_t* v)
{
  ASSERT(m > 0.0);

  body_t* b = polymec_malloc(sizeof(body_t));
  b->name = string_dup(name);
  b->m = m;
  b->x = *x;
  b->v = *v;
  return b;
}

void body_free(body_t* b)
{
  string_free(b->name);
  polymec_free(b);
}

