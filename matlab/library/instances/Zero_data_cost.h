#pragma once

class Zero_data_cost
{
public:
  Zero_data_cost(
    const matrix<double>& data,
    const matrix<int>& connectivity,
    const std::vector<double>& voxel_dimensions
  ) 
  {};

  template<typename R>
  R operator()(const R* const point1, const R* const point2) const
  {
    return R(0);
  }
};