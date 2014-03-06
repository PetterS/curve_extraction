#pragma once

template<typename R>
class Zero_data_cost
{
public:
  Zero_data_cost(
    const matrix<double>& data,
    const matrix<int>& connectivity,
    const std::vector<double>& voxel_dimensions
  ) 
  {};

  R operator () (R x1, R y1, R z1, R x2, R y2, R z2)
  {
    return R(0);
  }
};