#ifndef CURVE_SEGMENTATION_H
#define CURVE_SEGMENTATION_H

#include "mexutils.h"
#include "cppmatrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <tuple>
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


class length_cost_functor
{ 
  public:
    length_cost_functor (const vector<double>& vd_, float p_) 
      : voxeldimensions(vd_), penalty(p_) {};

    float operator () (float x1,float y1,float z1, float x2, float y2, float z2)
    {
      if (penalty == 0)
      {
        return 0;
      } else
      {
        float dx =voxeldimensions[0]*(x2-x1);
        float dy =voxeldimensions[1]*(y2-y1);
        float dz =voxeldimensions[2]*(z2-z1);

        return penalty*std::sqrt( dx*dx + dy*dy + dz*dz );
      }
    }

  private:
    vector<double> voxeldimensions;
    float penalty;
};

class curvature_cost_functor
{
  public: 
    curvature_cost_functor (const vector<double>& vd_, double p_, double pow_) 
      : voxeldimensions(vd_), penalty(p_), power(pow_) {};

  float operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3) 
  {
    return penalty* compute_curvature<float>
        (x1*voxeldimensions[0],y1*voxeldimensions[1],z1*voxeldimensions[2],
         x2*voxeldimensions[0],y2*voxeldimensions[1],z2*voxeldimensions[2],
         x3*voxeldimensions[0],y3*voxeldimensions[1],z3*voxeldimensions[2],
         power);
  }

  private:
    // Voxel dimensions
    const vector<double> voxeldimensions;
    double penalty;
    double power;
};

class torsion_cost_functor
{
  public: 
    torsion_cost_functor (const std::vector<double>& vd_, double p_, double pow_) 
      : voxeldimensions(vd_), penalty(p_), power(pow_) {};

  float operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3,
                      double x4, double y4, double z4) 
  {
    if (penalty == 0)
    {
      return 0;
    } else
    {
      return penalty* compute_torsion<float>(
          x1*voxeldimensions[0], y1*voxeldimensions[1], z1*voxeldimensions[2],
          x2*voxeldimensions[0], y2*voxeldimensions[1], z2*voxeldimensions[2],
          x3*voxeldimensions[0], y3*voxeldimensions[1], z3*voxeldimensions[2],
          x4*voxeldimensions[0], y4*voxeldimensions[1], z4*voxeldimensions[2],
          power);
    }
  }

  private:
    // Voxel dimensions
    const std::vector<double> voxeldimensions;
    double penalty;
    double power;
};


struct InstanceSettings
{
  InstanceSettings() :
         length_penalty(0),
         curvature_penalty(0),
         torsion_penalty(0),
         curvature_power(2.0),
         torsion_power(2.0),
         regularization_radius(4.0),
         use_a_star(false),
         verbose(false),
         unary_type("linear"),
         store_visit_time(false)
  { }

  double length_penalty;
  double curvature_penalty;
  double torsion_penalty;

  double curvature_power;
  double torsion_power;

  double regularization_radius;

  bool use_a_star;
  bool verbose;
  string unary_type;
  bool store_visit_time;
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
  settings.unary_type = params.get<string>("unary_type", "linear");

  if (settings.unary_type != "linear") 
    throw runtime_error("Only linear data term allowed");

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

float distance_between_points(float x1,float y1,float z1, float x2, float y2, float z2, const std::vector<double>& voxeldimensions);

void edge_segmentation( std::vector<Mesh::Point>& points,
                        double& run_time,
                        int& evaluations,
                        double& cost,
                        const matrix<int>& mesh_map,
                        PieceWiseConstant& data_term,
                        const matrix<int>& connectivity,
                        const InstanceSettings& settings,
                        const PointSets& start_sets,
                        const PointSets& end_sets,
                        const std::vector<double>& voxeldimensions,
                        const ShortestPathOptions& options,
                        matrix<double>& visit_time);

void  edgepair_segmentation( std::vector<Mesh::Point>& points,
                              double& run_time,
                              int& evaluations,
                              double& cost,
                              const matrix<int>& mesh_map,
                              PieceWiseConstant& data_term,
                              const matrix<int>& connectivity,
                              InstanceSettings& settings,
                              const std::vector<double>& voxeldimensions,
                              const ShortestPathOptions& options,
                              matrix<double>& visit_time
                             );
		 
void node_segmentation(std::vector<Mesh::Point>& points,
                      double& run_time,
                      int& evaluations,
                      double& cost,
                      const matrix<int>& mesh_map,
                      PieceWiseConstant& data_term,
                      const matrix<int>& connectivity,
                      const InstanceSettings& settings,
                      const PointSets& start_sets,
                      const PointSets& end_sets,
                      const std::vector<double>& voxeldimensions,
                      const ShortestPathOptions& options,
                      matrix<double>& visit_time
                      );

#endif
