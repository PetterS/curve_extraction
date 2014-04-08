#pragma once

class Euclidean_curvature
{
  public:
    Euclidean_curvature (
      const matrix<double>& data, 
      const InstanceSettings& settings)
      : dims(settings.voxel_dimensions), 
        penalty(settings.penalty[1]),
        power(settings.power[1]),
        data_dependent(false) 
  {};

  template<typename R>
  R operator()(const R* const point1, const R* const point2, const R* const point3) const
  {
    if (penalty == 0)
    {
       return R(0);
    } else
    {
    return  penalty * compute_curvature<R>(
            point1[0]*dims[0], point1[1]*dims[1], point1[2]*dims[2],
            point2[0]*dims[0], point2[1]*dims[1], point2[2]*dims[2],
            point3[0]*dims[0], point3[1]*dims[1], point3[2]*dims[2],
            power, false);
    }
  }

  const vector<double> dims;
  double penalty;
  double power;
  bool data_dependent;  
};