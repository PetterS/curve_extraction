class Linear_data_cost
{
  public: 
    Linear_data_cost(
          const matrix<double>& data, 
          const matrix<int>& connectivity,
          const std::vector<double>& voxel_dimensions) :
          data_term(data.data, data.M, data.N, data.O, voxel_dimensions)
  {};

  double operator () (double x1,double y1,double z1, double x2, double y2, double z2) 
  {
    return data_term.evaluate_line_integral<double>(x1,y1,z1, x2,y2,z2);
  }

protected:
  PieceWiseConstant data_term;
};