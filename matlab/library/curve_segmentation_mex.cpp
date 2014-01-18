// Johannes Ulén and Petter Strandmark 2013
#include "curve_segmentation.h"
bool verbose;
double timer;

#ifdef USE_OPENMP
#include <omp.h>
double get_wtime()
{
  return ::omp_get_wtime();
}
#else
#include <ctime>
double get_wtime()
{
  return std::time(0);
}
#endif

int M = 1;
int N = 1;
int O = 1;

// Work linear indices like MatLab, but starting from 0.
// Syntax coordinates (n1,n2,n3), image size (M,N,O);
bool validind(int n1, int n2, int n3)
{
  if ( (n1 > M-1 || n2 > N-1 || n3 > O-1) || (n1 < 0 || n2 < 0 || n3 < 0) )
    return false;

  return true;
}

bool validind(Mesh::Point p)
{
  return validind(p.x,p.y,p.z);
}

// Syntax coordinates (n1,n2,n3), image size (M,N,O);
int sub2ind(int n1, int n2, int n3)
{
  // Linear index
    return  n1 + n2*M + n3*M*N;
}

int sub2ind(Mesh::Point p)
{
  return sub2ind(p.x, p.y, p.z);
}

std::tuple<int,int,int> ind2sub(int n)
{
  int z = n/(M*N);
  int y = (n-z*M*N)/M;
  int x = n - y*M - z*M*N;

  return std::make_tuple(x,y,z);
}

Mesh::Point make_point(int n)
{
  int z = n/(M*N);
  int y = (n-z*M*N)/M;
  int x = n - y*M - z*M*N;

  return Mesh::Point(x,y,z);
}

void startTime()
{
  timer = ::get_wtime();
}

double endTime()
{
  double current_time = ::get_wtime();
  double elapsed = current_time - timer;
  timer = current_time;

  return elapsed;
}

double endTime(const char* message)
{
  double t = endTime();
  mexPrintf("%s : %g (s). \n", message, t);
  return t;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  startTime();

  // Check input and outputs
  ASSERT(nlhs == 6);
  ASSERT(nrhs == 3 || nrhs == 4); 
	
  // Mesh defining allowed pixels
  // Encoded as
  // 0: Disallowed
  // 1: Allowed
  // 2: Start set
  // 3: End set.
  int curarg =0;
  const matrix<unsigned char> mesh_map(prhs[curarg++]);
  const matrix<double> data(prhs[curarg++]);
  const matrix<int> connectivity(prhs[curarg++]);

  // For 2 images third column should be zeros.
  ASSERT(connectivity.N == 3);
  ASSERT(connectivity.ndim() == 2);

  M = mesh_map.M;
  N = mesh_map.N;
  O = mesh_map.O;

  // Only 2 or 3d grid
  if ((mesh_map.ndim() != 2) && (mesh_map.ndim() != 3))
      mexErrMsgTxt("Only two and three-dimensional problem supported. \n");

  // Optional options
  MexParams params(nrhs-curarg, prhs+curarg); //Structure to hold and parse additional parameters
  InstanceSettings settings = parse_settings(params); 

  // Check input
  ASSERT(settings.voxel_dimensions.size() == 3);
  ASSERT(settings.regularization_radius > 0);
  ASSERT(settings.length_penalty >= 0);
  ASSERT(settings.curvature_penalty >= 0);
  ASSERT(settings.torsion_penalty >= 0);

  // No torsion for 2D
  if ((mesh_map.ndim() == 2) && (settings.torsion_penalty != 0))
  {
    mexPrintf("Torsion is always zero in a plane. \n");
    settings.torsion_penalty = 0;
  }

  bool use_pairs = false;
  bool use_edges = false;

  // If torsion regularization is nonzero we need to use pairs
  if (settings.torsion_penalty != 0) {
    use_pairs = true;
  } else if (settings.curvature_penalty != 0)
  {
    use_edges = true;
  }

  if (verbose)
    mexPrintf("Connectivity size is %d. \n", connectivity.M);

  // Extra start and end sets. Cell arrays of points.
  const mxArray* start_sets_cell = params.get<const mxArray*>("start_sets", 0);
  const mxArray* end_sets_cell = params.get<const mxArray*>("end_sets", 0);
  PointSets start_sets, end_sets;
  if (start_sets_cell) {
    ASSERT(mxIsCell(start_sets_cell));
    auto cell_size = mxGetNumberOfElements(start_sets_cell);

    for (int i = 0; i < cell_size; ++i) {
      auto element = mxGetCell(start_sets_cell, i);
      matrix<double> points_matrix(element);

      start_sets.push_back(std::vector<Mesh::Point>());
      auto& points = start_sets.back();

      if (points_matrix.N == 3)
      {
        for (int j = 0; j < points_matrix.M; ++j) {
          points.push_back(Mesh::Point(points_matrix(j, 0),
                                       points_matrix(j, 1),
                                       points_matrix(j, 2)));
        }
      } else if (points_matrix.N == 2)
      {
        for (int j = 0; j < points_matrix.M; ++j) {
          points.push_back(Mesh::Point(points_matrix(j, 0),
                                       points_matrix(j, 1),
                                       0));
        }
      }
      else {
        mexErrMsgTxt("Error in defined start sets.");
      }
    }
  }
  if (end_sets_cell) {
    ASSERT(mxIsCell(end_sets_cell));
    auto cell_size = mxGetNumberOfElements(end_sets_cell);
    for (int i = 0; i < cell_size; ++i) {
      auto element = mxGetCell(end_sets_cell, i);
      end_sets.push_back(std::vector<Mesh::Point>());
      auto& points = end_sets.back();
      matrix<double> points_matrix(element);

      if (points_matrix.N == 3)
      {
        for (int j = 0; j < points_matrix.M; ++j) {
          points.push_back(Mesh::Point(points_matrix(j, 0),
                                       points_matrix(j, 1),
                                       points_matrix(j, 2)));
        }
      } else if (points_matrix.N == 2)
      {
        for (int j = 0; j < points_matrix.M; ++j) {
          points.push_back(Mesh::Point(points_matrix(j, 0),
                                       points_matrix(j, 1),
                                       0));
        }
      }
      else {
        mexErrMsgTxt("Error in defined end sets.");
      }
    }
  }

  if (verbose)
    endTime("Reading data");

  #ifdef USE_OPENMP
    int max_threads = omp_get_max_threads();
    if (settings.num_threads > 0) {
      max_threads = settings.num_threads;
    }
    omp_set_num_threads(max_threads);
    int current_num_threads = -1;
    if (verbose) {
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
  std::vector<int> empty_dimensions(3,0);
  std::vector<int> real_dimensions(3);
  std::vector<int> dimensions(3);

  real_dimensions[0] = mesh_map.M;
  real_dimensions[1] = mesh_map.N;
  real_dimensions[2] = mesh_map.O;

  if (options.store_visited)
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

  Data_cost data_cost(data, connectivity, settings);

  if (options.store_parents)
  {
  	// We want the entire tree.
  	options.compute_all_distances = true;
  }

  double run_time;
  double cost;
  int evaluations;
  std::vector<Mesh::Point> points;

  if (verbose)
  {
    mexPrintf("Regularization coefficients. Length: %g Curvature: %g Torsion: %g \n",
              settings.length_penalty, settings.curvature_penalty, settings.torsion_penalty);
    mexPrintf("Regularization powers: curvature: %g torsion %g \n",
              settings.curvature_power, settings.torsion_power);
  }

  // What kind of variables will be used in the graph?
  // Torsion: Pair of edges.
  // Curvature: Edges.
  // Length: Nodes.

  // Curvature and Length can be calculated on Pair of Edges but this is overkill.
  // Same goes for Length on edges.
  if (use_pairs)
  {
    edgepair_segmentation(points, run_time, evaluations, cost,
                          mesh_map, data_cost, connectivity, settings,
                          settings.voxel_dimensions, options, o_visit_map);
  }
  else if (use_edges)
  {
    edge_segmentation(  points, run_time, evaluations, cost,
                        mesh_map, data_cost, connectivity,
                        settings,
                        start_sets, end_sets, 
                        settings.voxel_dimensions, options, o_visit_map);
  } else 
  {
    node_segmentation( points, run_time, evaluations, cost,
                       mesh_map, data_cost,  connectivity,
                       settings, start_sets, end_sets,
                       settings.voxel_dimensions, options, 
                       o_visit_map, o_shortest_path_tree);
  } 

  matrix<double>  o_path(points.size(),3);
  matrix<double>  o_time(1);
  matrix<int>     o_eval(1);
  matrix<double>  o_cost(1);

  int n_line = 0;
  for (   auto it = points.begin(); it != points.end(); it++)
  {
      o_path(n_line, 0) = it->x;
      o_path(n_line, 1) = it->y;
      o_path(n_line, 2) = it->z;
      n_line++;
  }

  o_time(0) = run_time;
  o_eval(0) = evaluations;
  o_cost(0) = cost;

  // Write to MatLab
  plhs[0] = o_path;
  plhs[1] = o_cost;
  plhs[2] = o_time;
  plhs[3] = o_eval;
  plhs[4] = o_visit_map;
  plhs[5] = o_shortest_path_tree;  
}
