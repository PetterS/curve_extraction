#include <assert.h>
class Geodesic_length_cost
{
  public:
    Geodesic_length_cost (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty)
      : data(data),
        voxel_dimensions(voxel_dimensions), 
        penalty(penalty),
        data_depdent(true)
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

        assert(z1 == 0);
        assert(z2 == 0);

        // This is approximately true if |dx| <= 1, |dy| <= 1.
        double dz = voxel_dimensions[2]*(data(x2,y2)-data(x1,y1));
        return penalty*std::sqrt( dx*dx + dy*dy + dz*dz );
      }
    }

  bool data_depdent;     
  const matrix<double> data;
  vector<double> voxel_dimensions;
  double penalty;
};