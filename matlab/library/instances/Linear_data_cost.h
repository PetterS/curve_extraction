#pragma once

class Linear_data_cost
{
  public: 
    Linear_data_cost(
          const matrix<double>& data, 
          const matrix<int>& connectivity,
          const InstanceSettings& settings) :
          data_term(data.data, data.M, data.N, data.O, settings.voxel_dimensions)
  {};

  template<typename R>
  R operator()(const R* const point1, const R* const point2) const
  {
    return data_term.evaluate_line_integral(point1[0],point1[1],point1[2],
                                            point2[0],point2[1],point2[2]);
  }

protected:
  const PieceWiseConstant data_term;
};