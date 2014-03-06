#pragma once

template<typename R>
class Geodesic_curvature
{
  public:
    Geodesic_curvature (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty, 
      double power) :
    data_depdent(true)
  {};

  R operator () ( R x1, R y1, R z1,
                  R x2, R y2, R z2,
                  R x3, R y3, R z3)
  {
    return R(0);
  }

  bool data_depdent;  
};