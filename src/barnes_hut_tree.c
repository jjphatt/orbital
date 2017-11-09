// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#include "barnes_hut_tree.h"

struct barnes_hut_tree_t 
{
  MPI_Comm comm;
  bbox_t bbox;
  real_t theta;
  void (*compute_force)(void* context, 
                        int i, point_t* xi, 
                        int j, point_t* xj, 
                        vector_t* Fij);
};

barnes_hut_tree_t* barnes_hut_tree_new(MPI_Comm comm, 
                                       real_t theta,
                                       void (*compute_force)(void* context, int i, point_t* xi, int j, point_t* xj, vector_t* Fij))
{
  ASSERT(theta >= 0.0);
  ASSERT(compute_force != NULL);

  barnes_hut_tree_t* tree = polymec_malloc(sizeof(barnes_hut_tree_t));
  tree->comm = comm;
  tree->theta = theta;
  tree->compute_force = compute_force;
  return tree;
}

void barnes_hut_tree_free(barnes_hut_tree_t* tree)
{
  polymec_free(tree);
}

void barnes_hut_tree_compute_forces(barnes_hut_tree_t* tree,
                                    body_array_t* bodies,
                                    vector_t* forces)
{
  // Compute a bounding box that contains all of the points.
  bbox_t bbox;
  bbox_make_empty_set(&bbox);
  for (size_t i = 0; i < bodies->size; ++i)
    bbox_grow(&bbox, &(bodies->data[i]->x));

  // Insert the positions of the bodies into an octree.
  octree_t* octree = octree_new(&bbox);
  for (size_t i = 0; i < bodies->size; ++i)
    octree_insert(octree, &bodies->data[i]->x, (int)i);

  // Clean up.
  octree_free(octree);
}

