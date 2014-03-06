#pragma once

template<typename R>
class Linear_data_cost
{
  public: 
    Linear_data_cost(
          const matrix<double>& data, 
          const matrix<int>& connectivity,
          const std::vector<double>& voxel_dimensions) :
          data_term(data.data, data.M, data.N, data.O, voxel_dimensions)
  {};

  R operator () (R x1, R y1, R z1, R x2, R y2, R z2) 
  {
    return data_term.evaluate_line_integral<R>(x1,y1,z1, x2,y2,z2);
  }

protected:
  PieceWiseConstant data_term;
};