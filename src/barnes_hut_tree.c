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
  real_t total_mass;
  point_t center_of_mass;
  int num_points;
} agg_t;

typedef struct
{
  real_t G, theta;
  real_t Lx, Ly, Lz;
  agg_t* branches;
  real_t* masses;  // locally stored in the tree

  // Data for "test mass"
  real_t mi;
  point_t* xi;
  vector_t* Fi; 

  // Performance / diagnostics.
  int num_force_evals;
} force_calc_t;

// This function aggregates branch data to its parent branch.
static bool aggregate_branch(void* context, int depth, int branch_index, int parent_index)
{
  START_FUNCTION_TIMER();
  // Aggregate this branch to its parent.
  if (parent_index != -1)
  {
    force_calc_t* force_calc = context;
    agg_t* b_agg = &(force_calc->branches[branch_index]);
    agg_t* p_agg = &(force_calc->branches[parent_index]);
    p_agg->total_mass += b_agg->total_mass;
    p_agg->center_of_mass.x += b_agg->center_of_mass.x;
    p_agg->center_of_mass.y += b_agg->center_of_mass.y;
    p_agg->center_of_mass.z += b_agg->center_of_mass.z;
    p_agg->num_points += b_agg->num_points;
  }

  // Visit our children.
  STOP_FUNCTION_TIMER();
  return true;
}

// This function aggregates the children attached to a branch.
static void aggregate_leaf(void* context, int index, point_t* point, int parent_index)
{
  START_FUNCTION_TIMER();
  force_calc_t* force_calc = context;

  // Aggregate this point.
  agg_t* agg = &(force_calc->branches[parent_index]);
  agg->total_mass += force_calc->masses[index];
  agg->center_of_mass.x += point->x;
  agg->center_of_mass.y += point->y;
  agg->center_of_mass.z += point->z;
  ++(agg->num_points);
  STOP_FUNCTION_TIMER();
}

// This function accumulates forces from a branch.
static bool sum_branch_force(void* context, int depth, int branch_index, int parent_index)
{
  force_calc_t* force_calc = context;
  agg_t* agg = &(force_calc->branches[branch_index]);

  // If we have no points beneath us, we have nothing to do here.
  if (agg->num_points == 0)
    return false;

  // Divide our CM by the number of leaves beneath this branch.
  point_t* xj = &(agg->center_of_mass);
  xj->x /= agg->num_points;
  xj->y /= agg->num_points;
  xj->z /= agg->num_points;

  // Should this branch be aggregated?
  real_t L = MAX(force_calc->Lx, MAX(force_calc->Ly, force_calc->Lz));
  real_t s = L * pow(0.5, depth);
  point_t* xi = force_calc->xi;
  real_t d = point_distance(xi, xj);
  if ((s/d) < force_calc->theta)
  {
    // Yup! We compute the force as though this branch and its children are 
    // one big mass at the center of mass of the branch.
    real_t G = force_calc->G;
    real_t mi = force_calc->mi;
    real_t mj = agg->total_mass;
    vector_t rij;
    point_displacement(xi, xj, &rij);
    real_t r3_inv = 1.0 / (d*d*d);
    force_calc->Fi->x += G * mi * mj * rij.x * r3_inv;
    force_calc->Fi->y += G * mi * mj * rij.y * r3_inv;
    force_calc->Fi->z += G * mi * mj * rij.z * r3_inv;
    ++(force_calc->num_force_evals);
    return false; // Don't visit our children.
  }
  else
    return true; // Visit children.
}

// This function accumulates forces from a leaf.
static void sum_leaf_force(void* context, int j, point_t* xj, int parent_index)
{
  force_calc_t* force_calc = context;

  point_t* xi = force_calc->xi;
  real_t G = force_calc->G;
  real_t mi = force_calc->mi;
  real_t mj = force_calc->masses[j];
  vector_t rij;
  point_displacement(xi, xj, &rij);
  real_t r = vector_mag(&rij);
  if (r > 1e-12) // xi != xj
  {
    real_t r3_inv = 1.0 / (r*r*r);
    force_calc->Fi->x += G * mi * mj * rij.x * r3_inv;
    force_calc->Fi->y += G * mi * mj * rij.y * r3_inv;
    force_calc->Fi->z += G * mi * mj * rij.z * r3_inv;
  }
  ++(force_calc->num_force_evals);
}

void barnes_hut_tree_compute_forces(barnes_hut_tree_t* tree,
                                    real_t G,
                                    point_t* points,
                                    real_t* masses,
                                    int N,
                                    vector_t* forces)
{
  START_FUNCTION_TIMER();
  ASSERT(G > 0.0);

  // Compute a bounding box that contains all of the local points.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0, 
                 .z1 = 0.0, .z2 = 1.0};
  for (int i = 0; i < N; ++i)
    bbox_grow(&bbox, &(points[i]));

  // Insert the points into an octree.
  octree_t* octree = octree_new(&bbox);
  for (size_t i = 0; i < N; ++i)
    octree_insert(octree, &(points[i]), (int)i);

  // Set up a context to keep the books for our octree.
  force_calc_t force_calc;
  force_calc.Lx = bbox.x2 - bbox.x1;
  force_calc.Ly = bbox.y2 - bbox.y1;
  force_calc.Lz = bbox.z2 - bbox.z1;
  force_calc.G = G;
  force_calc.theta = tree->theta;
  force_calc.branches = octree_new_branch_array(octree, sizeof(agg_t));
  force_calc.masses = masses;
  force_calc.num_force_evals = 0;

  // Compute aggregates for points at each branch node of the octree.
  octree_visit(octree, OCTREE_PRE, &force_calc, 
               aggregate_branch, aggregate_leaf);

  // We compute forces for all points in our communicator, so let's get 
  // data from other processes.
  int nprocs, rank;
  MPI_Comm_size(tree->comm, &nprocs);
  MPI_Comm_rank(tree->comm, &rank);
  
  // How many points on other processes?
  int Np[nprocs];
  MPI_Allgather(&N, 1, MPI_INT, 
                Np, 1, MPI_INT, 
                tree->comm);
  int Ntot = 0;
  for (int p = 0; p < nprocs; ++p)
    Ntot += Np[p];

  // Get point and mass data from other processes.
  point_t* all_points = polymec_malloc(sizeof(point_t) * Ntot);
  MPI_Allgather(points, 3*N, MPI_REAL_T, 
                all_points, 3*N, MPI_REAL_T, 
                tree->comm);
  real_t* all_masses = polymec_malloc(sizeof(real_t) * Ntot);
  MPI_Allgather(masses, N, MPI_REAL_T, 
                all_masses, N, MPI_REAL_T, 
                tree->comm);

  // Now compute forces on all of the points. Note that we use a "post" 
  // traversal for this one, since we need to check whether the branches 
  // are aggregated.
  vector_t* all_forces = polymec_malloc(sizeof(vector_t) * Ntot);
  memset(all_forces, 0, sizeof(vector_t) * Ntot);
  for (int i = 0; i < Ntot; ++i)
  {
    force_calc.xi = &(all_points[i]);
    force_calc.mi = all_masses[i];
    force_calc.Fi = &(all_forces[i]);
    octree_visit(octree, OCTREE_POST, &force_calc, 
                 sum_branch_force, sum_leaf_force);
  }
  log_detail("barnes_hut_tree: theta = %g => %d force evaluations.", 
             tree->theta, force_calc.num_force_evals);

  // Sum together all forces on all points over all processes.
  MPI_Allreduce(MPI_IN_PLACE, all_forces, 3*Ntot, 
                MPI_REAL_T, MPI_SUM, tree->comm);

  // Extract force data for this process.
  int start = 0;
  for (int p = 0; p < rank; ++p)
    start += Np[p];
  memcpy(forces, &(all_forces[start]), sizeof(vector_t) * N);

  // Clean up.
  polymec_free(all_forces);
  polymec_free(all_masses);
  polymec_free(all_points);
  polymec_free(force_calc.branches);
  octree_free(octree);
  STOP_FUNCTION_TIMER();
}

