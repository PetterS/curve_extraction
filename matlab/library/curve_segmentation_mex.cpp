// Johannes Ulén and Petter Strandmark 2013
#include "curve_segmentation.h"

// Avoid explicit instantiation.
#include "node_segmentation.h"
#include "edge_segmentation.h"
#include "edgepair_segmentaion.h"

// Calls main_function
#include "instances/mex_wrapper_shortest_path.h"

// Data_cost: Any function of two points.
// Pair_cost any function of two points.
// Triplet_cost any function of three points.
// Quad_cost any function of four points. 
template<typename Data_cost, typename Pair_cost, typename Triplet_cost, typename Quad_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  startTime();

  // Check input and outputs
  ASSERT(nrhs == 4 || nrhs == 5);

  // Mesh defines allowed pixels/voxel encoded as
  // 0: Disallowed
  // 1: Allowed
  // 2: Start set
  // 3: End set.
  int curarg =1;
  const matrix<unsigned char> mesh_map(prhs[curarg++]);
  const matrix<double> data(prhs[curarg++]);
  const matrix<int> connectivity(prhs[curarg++]);

  // For 2 images third column should be zeros.
  ASSERT(connectivity.N == 3);
  ASSERT(connectivity.ndim() == 2);

  M = mesh_map.M;
  N = mesh_map.N;
  O = mesh_map.O;

  // Only 2d or 3d grid
  if ((mesh_map.ndim() != 2) && (mesh_map.ndim() != 3))
      mexErrMsgTxt("Only two and three-dimensional problem supported. \n");

  MexParams params(nrhs-curarg, prhs+curarg); //Structure to hold and parse additional parameters
  InstanceSettings settings = parse_settings(params);

  bool use_pairs = false;
  bool use_edges = false;

  if (settings.penalty[2] != 0) {
    use_pairs = true;
  } else if (settings.penalty[1] != 0)
  {
    use_edges = true;
  }

  // No line graph needed hence A* will not be used.
  if (!use_edges && !use_pairs)
    settings.use_a_star = false;

  if (settings.verbose)
    mexPrintf("Connectivity size is %d. \n", connectivity.M);

  if (settings.verbose)
    endTime("Reading data");

  #ifdef USE_OPENMP
    int max_threads = omp_get_max_threads();
    if (settings.num_threads > 0) {
      max_threads = settings.num_threads;
    }
    omp_set_num_threads(max_threads);
    int current_num_threads = -1;
    if (settings.verbose) {
      #pragma omp parallel for
      for (int i = 0; i < 1000; ++i)
      {
        current_num_threads = omp_get_num_threads();
      }
      mexPrintf("Using OpenMP with %d threads (maximum %d).\n",
                current_num_threads, max_threads);
    }
  #endif

  ShortestPathOptions options;
  options.print_progress = false;
  options.maximum_queue_size = 1000 * 1000 * 1000;
  options.store_visited = settings.store_visit_time;
  options.store_parents = settings.store_parents;

  // Empty matrices if visit order or parents are not calculated.
  // For edge and edge pair type graphs the visit map is needed
  // to resolve conflicts when shortest_path tree is to be calculated.
  std::vector<int> empty_dimensions(3,0);
  std::vector<int> real_dimensions(3);
  std::vector<int> dimensions(3);

  real_dimensions[0] = mesh_map.M;
  real_dimensions[1] = mesh_map.N;
  real_dimensions[2] = mesh_map.O;

  if (options.store_visited ||
      (use_edges && options.store_parents) ||
      (use_pairs && options.store_parents)
     )
    dimensions = real_dimensions;
  else
    dimensions = empty_dimensions;

  matrix<int>  o_visit_map   ( dimensions[0],
                               dimensions[1],
                               dimensions[2]);

  if (options.store_parents)
    dimensions = real_dimensions;
  else
    dimensions = empty_dimensions;

  matrix<int> o_shortest_path_tree( dimensions[0],
                                    dimensions[1],
                                    dimensions[2]);

  if (settings.store_distances)
    dimensions = real_dimensions;
  else
    dimensions = empty_dimensions;

  matrix<double> o_distances( dimensions[0],
                              dimensions[1],
                              dimensions[2]);

  if (settings.compute_all_distances)
  	options.compute_all_distances = true;

  double run_time;
  double cost;
  int evaluations;
  std::vector<Point> points;

  if (settings.verbose)
  {
    mexPrintf("Regularization coefficients Pair: %g Triplet: %g Quad: %g. \n",
              settings.penalty[0], settings.penalty[1], settings.penalty[2]);
    mexPrintf("Regularization powers Pair: %g Triplet: %g Quad: %g. \n",
              settings.power[0], settings.power[1], settings.power[2]);
  }

  // What kind of variables will be used in the graph?
  // Quad: Pair of edges.
  // Triplet: Edges.
  // Pair: Nodes.
  SegmentationOutput output(points, run_time, evaluations, cost, o_visit_map, o_shortest_path_tree, o_distances);

  // Triplet and Pair can be calculated on Pair of Edges but this is overkill.
  // Same goes for Pair on edges.
  if (use_pairs)
  {
    edgepair_segmentation<Data_cost, Pair_cost, Triplet_cost, Quad_cost>
    (data, mesh_map, connectivity, settings, options, output);
  }
  else if (use_edges)
  {
    edge_segmentation<Data_cost, Pair_cost, Triplet_cost>
    (data, mesh_map, connectivity, settings, options, output);
  }
  else
  {
    node_segmentation<Data_cost, Pair_cost>
    (data, mesh_map, connectivity, settings, options, output);
  }

  matrix<double>  o_time(1);
  matrix<int>     o_eval(1);
  matrix<double>  o_cost(1);

  int max_dim = 2;
  if (O > 1)
    max_dim = 3;
  
  matrix<double>  o_path(points.size(),max_dim);

  int n_line = 0;
  for ( auto p : points)
  {
    // +1 matlab index
    for (int dim = 0; dim < max_dim; dim++)
      o_path(n_line, dim) = p[dim] + 1;

    n_line++;
  }

  for (int i = 0; i < o_shortest_path_tree.numel(); i++)
    o_shortest_path_tree(i)++;
  

  o_time(0) = output.run_time;
  o_eval(0) = output.evaluations;
  o_cost(0) = output.cost;

  // Write to MatLab
  plhs[0] = o_path;
  plhs[1] = o_cost;
  plhs[2] = o_time;
  plhs[3] = o_eval;
  plhs[4] = o_visit_map;
  plhs[5] = o_shortest_path_tree;
  plhs[6] = o_distances;
}