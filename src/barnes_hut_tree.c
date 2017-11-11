// Copyright (c) 2017, Jeffrey N. Johnson
// All rights reserved.

#include "barnes_hut_tree.h"

struct barnes_hut_tree_t 
{
  MPI_Comm comm;
  real_t theta;
};

barnes_hut_tree_t* barnes_hut_tree_new(MPI_Comm comm, 
                                       real_t theta)
{
  ASSERT(theta >= 0.0);

  barnes_hut_tree_t* tree = polymec_malloc(sizeof(barnes_hut_tree_t));
  tree->comm = comm;
  tree->theta = theta;
  return tree;
}

void barnes_hut_tree_free(barnes_hut_tree_t* tree)
{
  polymec_free(tree);
}

// These structs store metadata for branches and nodes in the Barnes-Hut tree's 
// underlying octree.
typedef struct
{
  real_t total_charge;
  point_t centroid;
  int num_points;
} agg_t;

DEFINE_ARRAY(agg_array, agg_t)

typedef struct
{
  real_t K, theta;
  real_t Lx, Ly, Lz;
  agg_array_t* branches;
  real_t* charges;
  vector_t* forces;
  int i;
  point_t* x;
} force_calc_t;

// This function aggregates branch data to its parent branch.
static bool aggregate_branch(void* context, int depth, int branch_index, int parent_index)
{
  force_calc_t* force_calc = context;

  // Aggregate this branch to its parent.
  agg_t* b_agg = &(force_calc->branches->data[branch_index]);
  agg_t* p_agg = &(force_calc->branches->data[parent_index]);
  p_agg->total_charge += b_agg->total_charge;
  p_agg->centroid.x += b_agg->centroid.x;
  p_agg->centroid.y += b_agg->centroid.y;
  p_agg->centroid.z += b_agg->centroid.z;
  p_agg->num_points += b_agg->num_points;

  // Visit our children.
  return true;
}

// This function aggregates the children attached to a branch.
static void aggregate_leaf(void* context, int index, point_t* point, int parent_index)
{
  force_calc_t* force_calc = context;

  // Aggregate this point.
  agg_t* agg = &(force_calc->branches->data[parent_index]);
  agg->total_charge += force_calc->charges[index];
  agg->centroid.x += point->x;
  agg->centroid.y += point->y;
  agg->centroid.z += point->z;
  ++(agg->num_points);
}

// This function accumulates forces from a branch.
static bool sum_branch_force(void* context, int depth, int branch_index, int parent_index)
{
  force_calc_t* force_calc = context;

  // Normalize the centroid for this branch node.
  agg_t* agg = &(force_calc->branches->data[parent_index]);
  point_t* xj = &(agg->centroid);
  xj->x /= agg->num_points;
  xj->y /= agg->num_points;
  xj->z /= agg->num_points;

  // Should this branch be aggregated?
  real_t L = MAX(force_calc->Lx, MAX(force_calc->Ly, force_calc->Lz));
  real_t s = L * pow(0.5, depth);
  point_t* xi = force_calc->x;
  real_t d = point_distance(xi, xj);
  if (s / d < force_calc->theta)
  {
    // Yup! We compute the force as though this branch and its children are 
    // one big charge at the centroid.
    int i = force_calc->i;
    real_t K = force_calc->K;
    real_t qi = force_calc->charges[i];
    real_t qj = agg->total_charge;
    vector_t rij;
    point_displacement(xi, xj, &rij);
    real_t r3_inv = 1.0 / (d*d*d);
    force_calc->forces[i].x += K * qi * qj * rij.x * r3_inv;
    force_calc->forces[i].y += K * qi * qj * rij.y * r3_inv;
    force_calc->forces[i].z += K * qi * qj * rij.z * r3_inv;
    return false; // Don't visit our children.
  }
  else
    return true; // Visit children.
}

// This function accumulates forces from a leaf.
static void sum_leaf_force(void* context, int j, point_t* xj, int parent_index)
{
  force_calc_t* force_calc = context;

  int i = force_calc->i;
  if (i == j) return;

  real_t K = force_calc->K;
  real_t qi = force_calc->charges[i];
  real_t qj = force_calc->charges[j];
  point_t* xi = force_calc->x;
  vector_t rij;
  point_displacement(xi, xj, &rij);
  real_t r = vector_mag(&rij);
  real_t r3_inv = 1.0 / (r*r*r);
  force_calc->forces[i].x += K * qi * qj * rij.x * r3_inv;
  force_calc->forces[i].y += K * qi * qj * rij.y * r3_inv;
  force_calc->forces[i].z += K * qi * qj * rij.z * r3_inv;
}

void barnes_hut_tree_compute_forces(barnes_hut_tree_t* tree,
                                    real_t K,
                                    point_t* points,
                                    real_t* charges,
                                    int N,
                                    vector_t* forces)
{
  ASSERT(K > 0.0);

  // In the case of a single body, the force is zero.
  if (N == 1)
  {
    forces[0].x = forces[0].y = forces[0].z = 0.0;
    return;
  }

  // Compute a bounding box that contains all of the points.
  bbox_t bbox;
  bbox_make_empty_set(&bbox);
  for (int i = 0; i < N; ++i)
    bbox_grow(&bbox, &(points[i]));

  // Insert the points into an octree.
  octree_t* octree = octree_new(&bbox);
  for (size_t i = 0; i < N; ++i)
    octree_insert(octree, &(points[i]), (int)i);

  // Compute aggregates for points at each branch node of the octree.
  force_calc_t force_calc;
  force_calc.Lx = bbox.x2 - bbox.x1;
  force_calc.Ly = bbox.y2 - bbox.y1;
  force_calc.Lz = bbox.z2 - bbox.z1;
  force_calc.K = K;
  force_calc.theta = tree->theta;
  force_calc.branches = octree_new_branch_array(octree, sizeof(agg_t));
  force_calc.charges = charges;
  force_calc.forces = forces;
  octree_visit(octree, OCTREE_PRE, &force_calc, 
               aggregate_branch, aggregate_leaf);

  // Now compute forces on each of the points. Note that we use a "post" 
  // traversal for this one, since we need to check whether the branches 
  // are aggregated.
  for (int i = 0; i < N; ++i)
  {
    force_calc.i = i;
    force_calc.x = &(points[i]);
    octree_visit(octree, OCTREE_POST, &force_calc, 
                 sum_branch_force, sum_leaf_force);
  }

  // Clean up.
  polymec_free(force_calc.branches);
  octree_free(octree);
}

