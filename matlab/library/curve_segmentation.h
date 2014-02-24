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
typedef std::vector<std::vector<Mesh::Point>> PointSets;
extern double timer;
extern int M,N,O;
extern bool verbose;

enum Descent_method {lbfgs, nelder_mead};

struct InstanceSettings
{
  InstanceSettings()
  { }

  double length_penalty;
  double curvature_penalty;
  double torsion_penalty;

  double curvature_power;
  double torsion_power;

  double regularization_radius;

  bool use_a_star;
  bool verbose;

  bool store_visit_time;
  bool store_parents;
  bool store_distances;

  bool compute_all_distances;

  string data_type_str;

  vector<double> voxel_dimensions;

  double function_improvement_tolerance;
  double argument_improvement_tolerance;
  int num_threads;
  int maxiter;

  Descent_method descent_method;
  string descent_method_str;
};

InstanceSettings parse_settings(MexParams params)
{
  InstanceSettings settings;

  settings.verbose = params.get<bool>("verbose",false); // Debug messages
  settings.regularization_radius = params.get<double>("regularization_radius", 4.0);

  // Regularization coefficients
  settings.length_penalty    = params.get<double>("length_penalty", 0.0);
  settings.curvature_penalty = params.get<double>("curvature_penalty", 0.0);
  settings.torsion_penalty   = params.get<double>("torsion_penalty", 0.0);

  // Regularization is (curvature)^curvature_power.
  settings.curvature_power = params.get<double>("curvature_power", 2.0);
  settings.torsion_power = params.get<double>("torsion_power", 2.0);

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

  settings.descent_method_str = params.get<string>("descent_method","lbfgs");

  if (settings.descent_method_str == "lbfgs")
    settings.descent_method = lbfgs;
  else if (settings.descent_method_str == "nelder-mead")
    settings.descent_method = nelder_mead;
  else
    throw runtime_error("Unknown descent_method");

  settings.data_type_str = params.get<string>("data_type", "linear_interpolation");
  settings.voxel_dimensions = params.get< vector<double> >("voxel_dimensions");

  if (settings.voxel_dimensions.empty())
  {
    settings.voxel_dimensions.push_back(1.0);
    settings.voxel_dimensions.push_back(1.0);
    settings.voxel_dimensions.push_back(1.0);
  }

  return settings;
}

double get_wtime();
bool validind(int n1, int n2, int n3);
bool validind(Mesh::Point p);

int sub2ind(int n1, int n2, int n3);
int sub2ind(Mesh::Point p);

std::tuple<int,int,int> ind2sub(int n);
Mesh::Point make_point(int n);

std::tuple<int, int, int>
points_in_a_edgepair(int edgepair_num, const matrix<int>& connectivity);
std::vector<Mesh::Point>
edgepath_to_points(const std::vector<int>& path, const matrix<int>& connectivity);
double distance_between_points(double x1,double y1,double z1, double x2, double y2, double z2, const std::vector<double>& voxel_dimensions);

struct SegmentationOutput
{
  SegmentationOutput( std::vector<Mesh::Point>& points,
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

  std::vector<Mesh::Point>& points;
  double& run_time;
  int& evaluations;
  double& cost;
  matrix<int>& visit_time;
  matrix<int>& shortest_path_tree;
  matrix<double>& distances;
};

template<typename Data_cost, typename Length_cost>
void node_segmentation( const matrix<double>& data,
                        const matrix<unsigned char>& mesh_map,
                        const matrix<int>& connectivity,
                        InstanceSettings& settings,
                        ShortestPathOptions& options,
                        SegmentationOutput& output
                        );

template<typename Data_cost, typename Length_cost, typename Curvature_cost>
void edge_segmentation( const matrix<double>& data,
                        const matrix<unsigned char>& mesh_map,
                        const matrix<int>& connectivity,
                        InstanceSettings& settings,
                        ShortestPathOptions& options,
                        SegmentationOutput& output
                        );

template<typename Data_cost, typename Length_cost, typename Curvature_cost, typename Torsion_cost>
void  edgepair_segmentation(  const matrix<double>& data,
                              const matrix<unsigned char>& mesh_map,
                              const matrix<int>& connectivity,
                              InstanceSettings& settings,
                              ShortestPathOptions& options,
                              SegmentationOutput& output
                             );

#include "instances/instances.h"
#endif