#pragma once

class Zero_length
{
  public:
    Zero_length (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty)
      : data_depdent(false)
   {};

    template<typename R>
    R operator()(const R* const point1, const R* const point2) const
    {
      return R(0);
    }

    bool data_depdent;     // Can we cache this cost or not?
};