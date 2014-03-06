#pragma once

template<typename R>
class Euclidean_torsion
{
  public:
    Euclidean_torsion (
      const matrix<double>& data, 
      const std::vector<double>& voxel_dimensions, 
      double penalty, 
      double power)
      : voxel_dimensions(voxel_dimensions), 
        penalty(penalty), 
        power(power),
        data_depdent(false) {};

  R operator () ( R x1, R y1, R z1,
                  R x2, R y2, R z2,
                  R x3, R y3, R z3,
                  R x4, R y4, R z4)
  {
    if (penalty == 0)
    {
      return R(0);
    } else
    {
      return penalty* compute_torsion<R>(
          x1*voxel_dimensions[0], y1*voxel_dimensions[1], z1*voxel_dimensions[2],
          x2*voxel_dimensions[0], y2*voxel_dimensions[1], z2*voxel_dimensions[2],
          x3*voxel_dimensions[0], y3*voxel_dimensions[1], z3*voxel_dimensions[2],
          x4*voxel_dimensions[0], y4*voxel_dimensions[1], z4*voxel_dimensions[2],
          power, true);
    }
  }

  bool data_depdent;     // Can we cache this cost or not?
  const std::vector<double> voxel_dimensions;
  double penalty;
  double power;
};