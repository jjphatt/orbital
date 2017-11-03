// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#include "n_squared.h"

typedef struct
{
  real_t G;
  body_array_t* bodies;
} nsq_t;

static void nsq_init(void* context, real_t t)
{
//  nsq_t* nsq = context;
}

static real_t nsq_max_dt(void* context, real_t t, char* reason)
{
//  nsq_t* nsq = context;
  return REAL_MAX;
}

static real_t nsq_advance(void* context, real_t max_dt, real_t t)
{
//  nsq_t* nsq = context;
  return max_dt;
}

static void nsq_finalize(void* context, int step, real_t t)
{
}

static void nsq_add_probe(void* context, void* probe_context)
{
}

static void nsq_dtor(void* context)
{
  nsq_t* nsq = context;
  body_array_free(nsq->bodies);
  polymec_free(nsq);
}

model_t* n_squared_new(real_t G,
                       body_array_t* bodies)
{
  ASSERT(G >= 0.0);

  nsq_t* nsq = polymec_malloc(sizeof(nsq_t));
  nsq->G = G;
  nsq->bodies = bodies;
  model_vtable vtable = {.init = nsq_init,
                         .max_dt = nsq_max_dt,
                         .advance = nsq_advance,
                         .finalize = nsq_finalize,
                         .add_probe = nsq_add_probe,
                         .dtor = nsq_dtor};
  return model_new("N-squared", nsq, vtable, MODEL_SERIAL);
}

