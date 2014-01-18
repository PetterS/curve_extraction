#ifndef CURVE_SEGMENTATION_H
#define CURVE_SEGMENTATION_H

#include "mexutils.h"
#include "cppmatrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <map>
#include <memory>

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


class Length_cost
{
  public:
    Length_cost (const vector<double>& vd_, double p_)
      : voxel_dimensions(vd_), penalty(p_) {};

    double operator () (double x1,double y1,double z1, double x2, double y2, double z2)
    {
      if (penalty == 0)
      {
        return 0;
      } else
      {
        double dx = voxel_dimensions[0]*(x2-x1);
        double dy = voxel_dimensions[1]*(y2-y1);
        double dz = voxel_dimensions[2]*(z2-z1);

        return penalty*std::sqrt( dx*dx + dy*dy + dz*dz );
      }
    }

  private:
    vector<double> voxel_dimensions;
    double penalty;
};

class Curvature_cost
{
  public:
    Curvature_cost (const vector<double>& vd_, double p_, double pow_)
      : voxel_dimensions(vd_), penalty(p_), power(pow_) {};

  double operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3)
  {
    if (penalty == 0)
    {
       return 0;
    } else
    {
      return penalty* compute_curvature<double>
          (x1*voxel_dimensions[0],y1*voxel_dimensions[1],z1*voxel_dimensions[2],
           x2*voxel_dimensions[0],y2*voxel_dimensions[1],z2*voxel_dimensions[2],
           x3*voxel_dimensions[0],y3*voxel_dimensions[1],z3*voxel_dimensions[2],
           power, true);
    }
  }

  private:
    // Voxel dimensions
    const vector<double> voxel_dimensions;
    double penalty;
    double power;
};

class Torsion_cost
{
  public:
    Torsion_cost (const std::vector<double>& vd_, double p_, double pow_)
      : voxel_dimensions(vd_), penalty(p_), power(pow_) {};

  double operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3,
                      double x4, double y4, double z4)
  {
    if (penalty == 0)
    {
      return 0;
    } else
    {
      return penalty* compute_torsion<double>(
          x1*voxel_dimensions[0], y1*voxel_dimensions[1], z1*voxel_dimensions[2],
          x2*voxel_dimensions[0], y2*voxel_dimensions[1], z2*voxel_dimensions[2],
          x3*voxel_dimensions[0], y3*voxel_dimensions[1], z3*voxel_dimensions[2],
          x4*voxel_dimensions[0], y4*voxel_dimensions[1], z4*voxel_dimensions[2],
          power, true);
    }
  }

  private:
    // Voxel dimensions
    const std::vector<double> voxel_dimensions;
    double penalty;
    double power;
};

enum Descent_method {lbfgs, nelder_mead};
enum Unary_type {linear_interpolation, edge};

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

  Unary_type data_type;
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

  // Whether the time of visits should be stored.
  settings.store_visit_time = params.get<bool>("store_visit_time", false);
  settings.store_parents = params.get<bool>("store_parents", false);

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
  
  if (settings.data_type_str == "linear_interpolation")
    settings.data_type = linear_interpolation;
  else if (settings.data_type_str == "edge")
    settings.data_type = edge;
  else
    throw runtime_error("Unknown data type");

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

// Pure virtual class used as base for all Data costs
class Data_cost_base 
{
public:
  virtual double operator ()  (double x1,double y1,double z1, double x2, double y2, double z2) = 0;
};

class Linear_interpolation_data_cost : public Data_cost_base
{
  public: 
    Linear_interpolation_data_cost(
          const matrix<double>& data, 
          const std::vector<double>& voxel_dimensions) :
          data_term(data.data, data.M, data.N, data.O, voxel_dimensions)
  {};

  double operator () (double x1,double y1,double z1, double x2, double y2, double z2) 
  {
    return data_term.evaluate_line_integral<double>(x1,y1,z1, x2,y2,z2);
  }

  PieceWiseConstant data_term;
};

class Edge_data_cost : public Data_cost_base
{
public:
  Edge_data_cost(
    const matrix<double>& data,
    const matrix<int>& connectivity
  ) : data(data), connectivity(connectivity)
  {
    dims = data.ndim() -1;

    // Create fast table from (dx,dy,dz) to index in connectivity.
    for (int i = 0; i < connectivity.M; i++)
    {
      int dx,dy,dz;
      dx = connectivity(i,0);
      dy = connectivity(i,1);

      if (dims == 3)
        dz = connectivity(i,2);
      else
        dz = 0;

      lookup[ std::tuple<int, int, int>(dx,dy,dz) ] = i;
    }
  };

  double operator () (double x1,double y1,double z1, double x2, double y2, double z2)
  {
    int dx,dy,dz;

    dx = (int) x2 - x1;
    dy = (int) y2 - y1;

    if (dims == 3)
    {
      dz = (int) z2 - z1;
      int index = lookup[std::tuple<int, int, int>(dx,dy,dz)];
      return data(x1,y1,z1, index);
    }
    else
    {
      dz = 0;
      int index = lookup[std::tuple<int, int, int>(dx,dy,dz)];
      return data(x1,y1, index);
    }
  }

private:
  std::map<std::tuple<int, int, int>, int> lookup;
  const matrix<double> data;
  const matrix<int> connectivity;
  unsigned char dims;
};

// This wrapper class is introduced so that the data_cost can be
// passed as reference and initialized even though it's a pure virtual.
class Data_cost
{
public:
  ~Data_cost()
  {
    delete ptr;
  }

  Data_cost(matrix<double> data,
            matrix<int> connectivity,
            InstanceSettings settings)
  {
    if (settings.data_type == linear_interpolation)
      ptr = new Linear_interpolation_data_cost(data, settings.voxel_dimensions);
    else if (settings.data_type == edge)
      ptr = new Edge_data_cost(data, connectivity);
  }

  double operator ()  (double x1,double y1,double z1, double x2, double y2, double z2)
  {
    return ptr->operator() (x1,y1,z1, x2, y2, z2);
  }

private:
  Data_cost_base * ptr;
};


void edge_segmentation( std::vector<Mesh::Point>& points,
                        double& run_time,
                        int& evaluations,
                        double& cost,
                        const matrix<unsigned char>& mesh_map,
                        Data_cost& data_cost,
                        const matrix<int>& connectivity,
                        const InstanceSettings& settings,
                        const PointSets& start_sets,
                        const PointSets& end_sets,
                        const std::vector<double>& voxel_dimensions,
                        ShortestPathOptions& options,
                        matrix<int>& visit_time,
                        matrix<int>& shortest_path_tree);

void  edgepair_segmentation( std::vector<Mesh::Point>& points,
                              double& run_time,
                              int& evaluations,
                              double& cost,
                              const matrix<unsigned char>& mesh_map,
                              Data_cost& data_cost,
                              const matrix<int>& connectivity,
                              InstanceSettings& settings,
                              const std::vector<double>& voxel_dimensions,
                              ShortestPathOptions& options,
                              matrix<int>& visit_time,
                              matrix<int>& shortest_path_tree
                             );

void node_segmentation(std::vector<Mesh::Point>& points,
                      double& run_time,
                      int& evaluations,
                      double& cost,
                      const matrix<unsigned char>& mesh_map,
                      Data_cost& data_cost,
                      const matrix<int>& connectivity,
                      const InstanceSettings& settings,
                      const PointSets& start_sets,
                      const PointSets& end_sets,
                      const std::vector<double>& voxel_dimensions,
                      const ShortestPathOptions& options,
                      matrix<int>& visit_time,
                      matrix<int>& shortest_path_tree
                      );

#endif
