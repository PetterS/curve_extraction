#pragma once

template<typename R>
class Euclidean_length
{
  public:
    Euclidean_length (
      const matrix<R>& data, 
      const vector<R>& voxel_dimensions, 
      R penalty)
      : voxel_dimensions(voxel_dimensions), 
        penalty(penalty),
        data_depdent(false)
   {};

    R operator () (R x1, R y1, R z1, R x2, R y2, R z2)
    {
      if (penalty == 0)
      {
        return R(0);
      } else
      {
        R dx = voxel_dimensions[0]*(x2-x1);
        R dy = voxel_dimensions[1]*(y2-y1);
        R dz = voxel_dimensions[2]*(z2-z1);

        return penalty*std::sqrt( dx*dx + dy*dy + dz*dz );
      }
    }

    bool data_depdent;     // Can we cache this cost or not?
    vector<double> voxel_dimensions;
    double penalty;
};