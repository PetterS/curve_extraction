class Zero_curvature_cost
{
  public:
    Zero_curvature_cost (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty, 
      double power) :
    data_depdent(false)
  {};

  double operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3)
  {
    return 0;
  }

  bool data_depdent;  
};