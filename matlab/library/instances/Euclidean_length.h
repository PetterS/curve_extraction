#pragma once
class Euclidean_length
{
  public:
    Euclidean_length (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty)
      : dims(voxel_dimensions), 
        penalty(penalty),
        data_depdent(false)
   {};

    template<typename R>
    R operator()(const R* const point1, const R* const point2) const
    {
      using std::sqrt;
      R dx = dims[0]*(point1[0] - point2[0]);
      R dy = dims[1]*(point1[1] - point2[1]);
      R dz = dims[2]*(point1[2] - point2[2]);

      return  penalty*sqrt(dx*dx + dy*dy + dz*dz);
    }

    vector<double> dims;
    double penalty;
    bool data_depdent;  
};