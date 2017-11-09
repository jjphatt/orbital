// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#ifndef ORBITAL_BARNES_HUT_TREE_H
#define ORBITAL_BARNES_HUT_TREE_H

#include "body.h"

// A Barnes-Hut tree is a specialized, distributed octree that provides a 
// near/far-field approximation scheme for the n-body problem.
typedef struct barnes_hut_tree_t barnes_hut_tree_t;

// Creates a new Barnes-Hut tree on the given MPI communicator, containing the 
// given array of bodies. Arguments are:
//  comm          -- The MPI communicator on which the points are defined.
//  theta         -- A parameter that determines when to use the far-field 
//                   approximation for a force between two bodies 
//                   (typically 0.5).
//  compute_force -- A function that computes the force on point i due to 
//                   point j, storing the result in the vector Fij.
// typically 0.5. The tree copies data from the array of bodies.
barnes_hut_tree_t* barnes_hut_tree_new(MPI_Comm comm, 
                                       real_t theta,
                                       void (*compute_force)(void* context, int i, point_t* xi, int j, point_t* xj, vector_t* Fij));

// Frees the Barnes-Hut tree.
void barnes_hut_tree_free(barnes_hut_tree_t* tree);

// Computes the forces on the given set of bodies, storing the force 
// vectors in the forces array. forces[i] contains the net force on the ith 
// body in bodies.
void barnes_hut_tree_compute_forces(barnes_hut_tree_t* tree,
                                    body_array_t* bodies,
                                    vector_t* forces);                                     

#endif

