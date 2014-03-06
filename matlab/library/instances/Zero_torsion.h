#pragma once

template<typename R>
class Zero_torsion
{
  public:
    Zero_torsion (
      const matrix<double>& data, 
      const std::vector<double>& voxel_dimensions, 
      double penalty, 
      double power)
      :  data_depdent(false) {};

  R operator () ( R x1, R y1, R z1,
                  R x2, R y2, R z2,
                  R x3, R y3, R z3,
                  R x4, R y4, R z4)
  {
    return R(0);
  }

  bool data_depdent;  
};