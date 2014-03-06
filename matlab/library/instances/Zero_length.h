#pragma once

template<typename R>
class Zero_length
{
  public:
    Zero_length (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty)
      : data_depdent(false)
   {};

    R operator () (R x1, R y1, R z1, R x2, R y2, R z2)
    {
      return R(0);
    }

    bool data_depdent;     // Can we cache this cost or not?
};