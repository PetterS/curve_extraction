// Johannes Ulén and Petter Strandmark 2013
#include "curve_segmentation.h"

// Avoid explicit instantiation.
#include "node_segmentation.cpp"
#include "edge_segmentation.cpp"
#include "edgepair_segmentaion.cpp"

bool verbose;
double timer;

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

// Data_cost: Any function of two points.
// Length_cost any function of two points.
// Curvature_cost any function of three points.
// Torsion_ocst any function of four points.
template<typename Data_cost, typename Length_cost, typename Curvature_cost, typename Torsion_cost>
void curve_segmentation(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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

  // No line graph needed hence A* will not be used.
  if (!use_edges && !use_pairs)
    settings.use_a_star = false;

  if (verbose)
    mexPrintf("Connectivity size is %d. \n", connectivity.M);

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
  SegmentationOutput output(points, run_time, evaluations, cost, o_visit_map, o_shortest_path_tree, o_distances);

  // Curvature and Length can be calculated on Pair of Edges but this is overkill.
  // Same goes for Length on edges.
  if (use_pairs)
  {
    edgepair_segmentation<Data_cost, Length_cost, Curvature_cost, Torsion_cost>
    (data, mesh_map, connectivity, settings, options, output);
  }
  else if (use_edges)
  {
    edge_segmentation<Data_cost, Length_cost, Curvature_cost>
    (data, mesh_map, connectivity, settings, options, output);
  }
  else
  {
    node_segmentation<Data_cost, Length_cost>
    (data, mesh_map, connectivity, settings, options, output);
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

// Wrapper data from MATLAB.
void mexFunction(int            nlhs,     /* number of expected outputs */
                 mxArray        *plhs[],  /* mxArray output pointer array */
                 int            nrhs,     /* number of inputs */
                 const mxArray  *prhs[]   /* mxArray input pointer array */)
{
 char problem_type[1024];
 if (mxGetString(prhs[0], problem_type, 1024))
   throw runtime_error("First argument must be a string.");

  if (!strcmp(problem_type,"linear_interpolation"))
    curve_segmentation<Linear_data_cost, Euclidean_length, Euclidean_curvature, Euclidean_torsion>(nlhs, plhs, nrhs, prhs);
  else if (!strcmp(problem_type,"edge"))
    curve_segmentation<Edge_data_cost, Euclidean_length, Euclidean_curvature, Euclidean_torsion>(nlhs, plhs, nrhs, prhs); 
  else if (!strcmp(problem_type,"geodesic"))
    curve_segmentation<Zero_data_cost, Geodesic_length, Geodesic_curvature, Zero_torsion>(nlhs, plhs, nrhs, prhs);
  else
    throw runtime_error("Unknown data type");
}