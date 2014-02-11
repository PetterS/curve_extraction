class Zero_data_cost 
{
public:
  Zero_data_cost(
    const matrix<double>& data,
    const matrix<int>& connectivity,
    const std::vector<double>& voxel_dimensions
  ) 
  {};

  double operator () (double x1,double y1,double z1, double x2, double y2, double z2)
  {
    return 0;
  }
};