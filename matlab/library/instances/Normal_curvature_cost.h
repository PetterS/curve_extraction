class Normal_curvature_cost
{
  public:
    Normal_curvature_cost (
      const matrix<double>& data, 
      const vector<double>& vd_, 
      double p_, 
      double pow_)
      : voxel_dimensions(vd_), penalty(p_), power(pow_) {};

  double operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3)
  {
    if (penalty == 0)
    {
       return 0;
    } else
    {
      return penalty* compute_curvature<double>
          (x1*voxel_dimensions[0],y1*voxel_dimensions[1],z1*voxel_dimensions[2],
           x2*voxel_dimensions[0],y2*voxel_dimensions[1],z2*voxel_dimensions[2],
           x3*voxel_dimensions[0],y3*voxel_dimensions[1],z3*voxel_dimensions[2],
           power, true);
    }
  }

  protected:
    const vector<double> voxel_dimensions;
    double penalty;
    double power;
};