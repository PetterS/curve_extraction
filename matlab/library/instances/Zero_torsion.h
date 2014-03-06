#pragma once

class Zero_torsion
{
  public:
    Zero_torsion (
      const matrix<double>& data, 
      const std::vector<double>& voxel_dimensions, 
      double penalty, 
      double power)
      :  data_depdent(false) {};

  template<typename R>    
  R operator()(const R* const point1, const R* const point2, const R* const point3, const R* const point4) const
  {
    return R(0);
  }

  bool data_depdent;  
};