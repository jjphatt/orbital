// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#ifndef ORBITAL_BARNES_HUT_TREE_H
#define ORBITAL_BARNES_HUT_TREE_H

#include "body.h"

// A Barnes-Hut tree is a specialized, distributed octree that provides a 
// near/far-field approximation scheme for the n-body problem. Specifically, 
// this tree computes the forces on a set of points, where the forces are of 
// the form
//
//                     K * qi * qj
//           F(i, j) = -----------
//                               2
//                      |xi - xj|
//
typedef struct barnes_hut_tree_t barnes_hut_tree_t;

// Creates a new Barnes-Hut tree on the given MPI communicator, containing the 
// given array of bodies. Arguments are:
//  comm          -- The MPI communicator on which the points are defined.
//  theta         -- A parameter that determines when to use the far-field 
//                   approximation for a force between two bodies 
//                   (typically 0.5).
barnes_hut_tree_t* barnes_hut_tree_new(MPI_Comm comm, 
                                       real_t theta);

// Frees the Barnes-Hut tree.
void barnes_hut_tree_free(barnes_hut_tree_t* tree);

// Computes the forces on the given N points, storing the force 
// vectors in the forces array. Arguments:
//  K       -- The coupling constant for the force F(i, j).
//  points  -- An array of N points.
//  charges -- An array of N charges for the points.
//  N       -- The number of points.
//  forces  -- An array with room to store N vector-valued forces.
void barnes_hut_tree_compute_forces(barnes_hut_tree_t* tree,
                                    real_t K,
                                    point_t* points,
                                    real_t* charges,
                                    int N,
                                    vector_t* forces);

#endif

