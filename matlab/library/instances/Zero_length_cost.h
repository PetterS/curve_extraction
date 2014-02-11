class Zero_length_cost
{
  public:
    Zero_length_cost (
      const matrix<double>& data, 
      const vector<double>& voxel_dimensions, 
      double penalty)
      : data_depdent(false)
   {};

    double operator () (double x1,double y1,double z1, double x2, double y2, double z2)
    {
      return 0;
    }

    bool data_depdent;     // Can we cache this cost or not?
};