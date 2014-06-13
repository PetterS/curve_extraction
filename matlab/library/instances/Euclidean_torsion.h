#pragma once

class Euclidean_torsion
{
  public:
    Euclidean_torsion (
      const matrix<double>& data, 
      const InstanceSettings& settings)
      : dims(settings.voxel_dimensions), 
        penalty(settings.penalty[2]),
        power(settings.power[2]),
        data_dependent(false) {};

  template<typename R>      
  R operator()(const R* const point1,
               const R* const point2,
               const R* const point3,
               const R* const point4) const
  {
    return  penalty * compute_torsion<R>(
            point1[0]*dims[0], point1[1]*dims[1], point1[2]*dims[2],
            point2[0]*dims[0], point2[1]*dims[1], point2[2]*dims[2],
            point3[0]*dims[0], point3[1]*dims[1], point3[2]*dims[2],
            point4[0]*dims[0], point4[1]*dims[1], point4[2]*dims[2],
            power, false);
  }

  const std::vector<double> dims;
  double penalty;
  double power;
  bool data_dependent;
};