#ifndef CURVE_SEGMENTATION_H
#define CURVE_SEGMENTATION_H

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

#include "mexutils.h"
#include "cppmatrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <map>
#include <memory>
#include <string.h>

using std::ignore;
using std::tie;

#include <curve_extraction/curvature.h>
#include <curve_extraction/data_term.h>
#include <curve_extraction/grid_mesh.h>
#include <curve_extraction/shortest_path.h>

using namespace curve_extraction;
double timer;
int M = 1;
int N = 1;
int O = 1;
const int max_index = std::numeric_limits<int>::max();
enum Descent_method {lbfgs, nelder_mead};

struct Point
{
  Point () {}
  Point (double x, double y, double z)
  {
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }

  bool operator== (Point &other)
  {
    return (xyz[0] == other.xyz[0] &&
            xyz[1] == other.xyz[1] &&
            xyz[2] == other.xyz[2]);
  }

  double& operator [](int i) {
    return xyz[i];
  }

  double xyz[3];
};

class Delta_point
{
  public:
    Delta_point (const matrix<int>& connectivity, bool reverse_connectivity = false) :
    connectivity(connectivity), reverse_connectivity(reverse_connectivity)
  {}

  Point operator () (Point& root, int k) 
  {
    Point p;
    if (!reverse_connectivity)
    {
      p[0] = root[0] + connectivity(k,0);
      p[1] = root[1] + connectivity(k,1);
      p[2] = root[2] + connectivity(k,2);
    } else
    {
      p[0] = root[0] - connectivity(k,0);
      p[1] = root[1] - connectivity(k,1);
      p[2] = root[2] - connectivity(k,2);   
    }

    return p;
  }

  Point reverse(Point& root, int k) 
  {
    Point p;
    if (!reverse_connectivity)
    {
      p[0] = root[0] - connectivity(k,0);
      p[1] = root[1] - connectivity(k,1);
      p[2] = root[2] - connectivity(k,2);
    } else
    {
      p[0] = root[0] + connectivity(k,0);
      p[1] = root[1] + connectivity(k,1);
      p[2] = root[2] + connectivity(k,2);   
    }
    return p;
  }

  int size()
  {
    return connectivity.M;
  }


protected:
  bool reverse_connectivity;
  const matrix<int> connectivity;  
};

struct InstanceSettings
{
  InstanceSettings() :
  penalty(4,0.0),
  power(4,1.0),
  local_limit(4, std::numeric_limits<double>::infinity())
  { }

  // penalty[0]: e.g. length
  // penalty[1]: e.g. curvature
  // penalty[2]: e.g. torsion
  std::vector<double> penalty;
  std::vector<double> power;
  std::vector<double> local_limit;

  bool use_a_star;
  bool verbose;

  bool store_visit_time;
  bool store_parents;
  bool store_distances;

  bool compute_all_distances;

  bool fully_contained_set;
  
  string data_type_str;

  vector<double> voxel_dimensions;

  double function_improvement_tolerance;
  double argument_improvement_tolerance;

  double length_local_limit;
  double curvature_local_limit;
  double torsion_local_limit;
  double jounce_local_limit;

  int num_threads;
  int maxiter;

  Descent_method descent_method;
  string descent_method_str;
};

InstanceSettings parse_settings(MexParams params)
{
  InstanceSettings settings;


  settings.verbose = params.get<bool>("verbose",false); // Debug messages

  settings.penalty = params.get< vector<double> >("penalty");
  settings.power = params.get< vector<double> >("power");
  settings.local_limit = params.get< vector<double> >("local_limit");
  settings.voxel_dimensions = params.get< vector<double> >("voxel_dimensions");

  for (double& p : settings.penalty)
    ASSERT(p >= 0);

  for (double& p : settings.power)
    ASSERT(p >= 0);

  for (double& p : settings.local_limit)
    ASSERT(p >= 0);

  // Whether A* should be used for curvature.
  settings.use_a_star = params.get<bool>("use_a_star", false);

  // Store visit time for each node.
  settings.store_visit_time = params.get<bool>("store_visit_time", false);

  // Store the parent to each node.
  settings.store_parents = params.get<bool>("store_parents", false);

  // Store distance to each node.
  settings.store_distances = params.get<bool>("store_distances", false);

  // Visit the full graph
  settings.compute_all_distances = params.get<bool>("compute_all_distances", false);

  // Used by local optimization
  settings.function_improvement_tolerance = params.get<double>("function_improvement_tolerance", 1e-12);
  settings.argument_improvement_tolerance = params.get<double>("argument_improvement_tolerance", 1e-12);
  settings.num_threads = params.get<int>("num_threads", -1);
  settings.maxiter = params.get<int>("maxiter", 1000);

  // Only add edges _fully_ contained in the start and end set
  // At the moment only used for dubins path
  settings.fully_contained_set = params.get<bool>("fully_contained_set", false);


  settings.descent_method_str = params.get<string>("descent_method","lbfgs");

  if (settings.descent_method_str == "lbfgs")
    settings.descent_method = lbfgs;
  else if (settings.descent_method_str == "nelder-mead")
    settings.descent_method = nelder_mead;
  else
    throw runtime_error("Unknown descent method");

  settings.data_type_str = params.get<string>("data_type", "linear_interpolation");


  return settings;
}

// Work linear indices like MatLab, but starting from 0.
// Syntax coordinates (n1,n2,n3), image size (M,N,O);
bool validind(int n1, int n2, int n3)
{
  if ( (n1 > M-1 || n2 > N-1 || n3 > O-1) || (n1 < 0 || n2 < 0 || n3 < 0) )
    return false;

  return true;
}

bool valid_point(Point p)
{
  return validind(p[0],p[1],p[2]);
}

// Syntax coordinates (n1,n2,n3), image size (M,N,O);
int sub2ind(int n1, int n2, int n3)
{
  // Linear index
    return  n1 + n2*M + n3*M*N;
}

int point2ind(Point p)
{
  return sub2ind(p[0], p[1], p[2]);
}

std::tuple<int,int,int> ind2sub(int n)
{
  int z = n/(M*N);
  int y = (n-z*M*N)/M;
  int x = n - y*M - z*M*N;

  return std::make_tuple(x,y,z);
}

Point make_point(int n)
{
  int z = n/(M*N);
  int y = (n-z*M*N)/M;
  int x = n - y*M - z*M*N;

  return Point(x,y,z);
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

struct SegmentationOutput
{
  SegmentationOutput( std::vector<Point>& points,
                      double& run_time,
                      int& evaluations,
                      double& cost,
                      matrix<int>& visit_time,
                      matrix<int>& shortest_path_tree,
                      matrix<double>& distances) :
    points(points), run_time(run_time), evaluations(evaluations),
    cost(cost), visit_time(visit_time), shortest_path_tree(shortest_path_tree),
    distances(distances)
  {};

  std::vector<Point>& points;
  double& run_time;
  int& evaluations;
  double& cost;
  matrix<int>& visit_time;
  matrix<int>& shortest_path_tree;
  matrix<double>& distances;
};

#include "instances/instances.h"
#endif