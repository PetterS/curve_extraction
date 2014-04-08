#pragma once

class Zero_pair
{
  public:
    Zero_pair (
      const matrix<double>& data, 
      const InstanceSettings& settings)
      : data_dependent(false)
   {};

    template<typename R>
    R operator()(const R* const point1, const R* const point2) const
    {
      return R(0);
    }

    bool data_dependent;     
};