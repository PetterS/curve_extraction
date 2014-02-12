class Euclidean_length
{
  public:
    Euclidean_length (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty)
      : voxel_dimensions(voxel_dimensions), 
        penalty(penalty),
        data_depdent(false)
   {};

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

    bool data_depdent;     // Can we cache this cost or not?
    vector<double> voxel_dimensions;
    double penalty;
};