// Copyright (c) 2016, Jeffrey N. Johnson
// All rights reserved.

#include "body.h"

static void body_free(void* context)
{
  body_t* b = context;
  string_free(b->name);
}

body_t* body_new(const char* name,
                 real_t m,
                 point_t* x, 
                 vector_t* v)
{
  ASSERT(m > 0.0);

  body_t* b = polymec_refcounted_malloc(sizeof(body_t), body_free);
  b->name = string_dup(name);
  b->m = m;
  b->x = *x;
  b->v = *v;
  return b;
}

body_t* body_clone(body_t* b)
{
  return body_new(b->name, b->m, &(b->x), &(b->v)); 
}

