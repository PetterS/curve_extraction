#pragma once
#include "Boundary_points.h"

typedef std::tuple<double,double,double,double> corners;

class Depth_corners
{
public:
  Depth_corners(const matrix<double>& data) : data(data)
  {}

  corners operator() (boundary_point bp)
  {
    return operator()(bp.first, bp.second);
  }

  corners operator () (double x,double y)
  {
    // Determine four closest points
    double x_int, y_int;
    modf(x, &x_int);
    modf(y, &y_int);

    double i00,i10,i01,i11;

    // Sample data_cost
    i00 = data_value(x_int+0,y_int+0);
    i01 = data_value(x_int+0,y_int+1);
    i10 = data_value(x_int+1,y_int+0);
    i11 = data_value(x_int+1,y_int+1);
     
    return corners(i00,i01,i10,11);
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
};