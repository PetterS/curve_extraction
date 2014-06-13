#pragma once

class Zero_pentuple
{
  public:
    Zero_pentuple (
      const matrix<double>& data, 
      const InstanceSettings& settings) :  
        data_dependent(false)
      {};

  template<typename R>    
  R operator()(const R* const point1, const R* const point2, const R* const point3, const R* const point4, const R* const point5) const
  {
    return R(0);
  }

  bool data_dependent;
};