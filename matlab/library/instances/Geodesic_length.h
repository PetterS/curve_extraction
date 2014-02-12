#pragma once
#include <assert.h>
#include "Boundary_points.h"
#include "Depth_function.h"

class Geodesic_length
{
  public:
    Geodesic_length (
      const matrix<double>& data,
      const vector<double>& voxel_dimensions,
      double penalty)
      : data(data),
        voxel_dimensions(voxel_dimensions),
        penalty(penalty),
        data_depdent(true),
        boundary_points(voxel_dimensions, data.M, data.N),
        euclidean_length(data, voxel_dimensions, penalty),
        depth(data, voxel_dimensions)
   {
   }

    double operator () (double x1,double y1,double z1, double x2, double y2, double z2)
    {
      double dx = voxel_dimensions[0]*(x2-x1);
      double dy = voxel_dimensions[1]*(y2-y1);

      // More costly interpolation for longer interactions
      // For fixed connectivity relative boundary_points can be cached.
      std::vector<boundary_point> points = boundary_points(x1,y1,x2,y2);
      double cost = 0;

      auto prev = points.begin();
      auto next = prev;
      next++;

      // Length on each plane of the surface.
      z1 = depth(*prev);

      for (; next != points.end(); prev++,next++)
      {
          z2 = depth(*next);
          cost += euclidean_length(prev->first,prev->second,z1,
                                   next->first,next->second,z2);
          z1 = z2;
      }

      return penalty*cost;
    }

  // Data
  bool data_depdent;
  const matrix<double> data;
  vector<double> voxel_dimensions;
  double penalty;

  // Functors
  Boundary_points boundary_points;
  Euclidean_length euclidean_length;
  Depth_function depth;
};