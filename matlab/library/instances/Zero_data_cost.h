#pragma once

class Zero_data_cost
{
public:
  Zero_data_cost(
    const matrix<double>& data,
    const matrix<int>& connectivity,
    const InstanceSettings& settings
  ) 
  {};

  template<typename R>
  R operator()(const R* const point1, const R* const point2) const
  {
    return R(0);
  }
};