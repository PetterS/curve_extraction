#pragma once
#include "Boundary_points.h"

class Depth_function
{
public:
  Depth_function(const matrix<double>& data, const vector<double>& voxel_dimensions) :
    data(data), vd(voxel_dimensions)
  {}

  double operator() (boundary_point bp)
  {
    return operator()(bp.first, bp.second);
  }

  double operator () (double x,double y)
  {
    // Determine four closest points
    double x_int, y_int;
    modf(x, &x_int);
    modf(y, &y_int);

    double i00,i10,i01,i11;

    // Sample data_cost
    i00 = data_value(x_int+0,y_int+0);
    i10 = data_value(x_int+1,y_int+0);
    i01 = data_value(x_int+0,y_int+1);
    i11 = data_value(x_int+1,y_int+1);
     
    // Location in the square
    x = x-x_int;
    y = y-y_int;

    // The square is defined by [0,vd[0]]x[0,vd[1]
    // Bilinear interpolation
    return  ( 1/ ((vd[0]-0)*(vd[1]-0) ) ) *
            (  i00 * (vd[0]-x) *(vd[1]-y)
             + i10 * (x-0)     *(vd[1]-y)
            + i01 * (vd[0]-x) *(y-0)
             + i11 * (x-0)     *(y-0)
            );
  } 

  // Sample data or closest point inside data.
  double data_value(int x, int y)
  {
    x = feasible(x, data.M);
    y = feasible(y, data.N);

    return data(x,y);
  }

  int feasible(int x, int M)
  {
    if (x < 0)
      return 0;

    if (x > M-1)
      return M-1;

    return x;
  }

  const matrix<double> data;
  const vector<double> vd;
};