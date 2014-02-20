#pragma once
class Zero_torsion
{
  public:
    Zero_torsion (
      const matrix<double>& data, 
      const std::vector<double>& voxel_dimensions, 
      double penalty, 
      double power)
      :  data_depdent(false) {};

  double operator () (double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3,
                      double x4, double y4, double z4)
  {
    return 0;
  }

  bool data_depdent;  
};