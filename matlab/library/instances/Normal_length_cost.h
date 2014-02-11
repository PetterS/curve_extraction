class Normal_length_cost
{
  public:
    Normal_length_cost (
      const matrix<double>& data, 
      const vector<double>& vd_, 
      double p_)
      : voxel_dimensions(vd_), penalty(p_) {};

    double operator () (double x1,double y1,double z1, double x2, double y2, double z2)
    {
      if (penalty == 0)
      {
        return 0;
      } else
      {
        double dx = voxel_dimensions[0]*(x2-x1);
        double dy = voxel_dimensions[1]*(y2-y1);
        double dz = voxel_dimensions[2]*(z2-z1);

        return penalty*std::sqrt( dx*dx + dy*dy + dz*dz );
      }
    }

  protected:
    vector<double> voxel_dimensions;
    double penalty;
};