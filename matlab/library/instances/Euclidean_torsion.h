class Euclidean_torsion
{
  public:
    Euclidean_torsion (
      const matrix<double>& data, 
      const std::vector<double>& voxel_dimensions, 
      double penalty, 
      double power)
      : voxel_dimensions(voxel_dimensions), 
        penalty(penalty), 
        power(power),
        data_depdent(false) {};

  double operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3,
                      double x4, double y4, double z4)
  {
    if (penalty == 0)
    {
      return 0;
    } else
    {
      return penalty* compute_torsion<double>(
          x1*voxel_dimensions[0], y1*voxel_dimensions[1], z1*voxel_dimensions[2],
          x2*voxel_dimensions[0], y2*voxel_dimensions[1], z2*voxel_dimensions[2],
          x3*voxel_dimensions[0], y3*voxel_dimensions[1], z3*voxel_dimensions[2],
          x4*voxel_dimensions[0], y4*voxel_dimensions[1], z4*voxel_dimensions[2],
          power, true);
    }
  }

  bool data_depdent;     // Can we cache this cost or not?
  const std::vector<double> voxel_dimensions;
  double penalty;
  double power;
};