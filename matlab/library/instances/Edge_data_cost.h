#pragma once
// Explicitly defined cost for each edge.
class Edge_data_cost 
{
public:
  Edge_data_cost(
    const matrix<double>& data,
    const matrix<int>& connectivity,
    const std::vector<double>& voxel_dimensions
  ) : data(data), connectivity(connectivity)
  {
    dims = data.ndim() -1;

    // Create fast table from (dx,dy,dz) to index in connectivity.
    for (int i = 0; i < connectivity.M; i++)
    {
      int dx,dy,dz;
      dx = connectivity(i,0);
      dy = connectivity(i,1);

      if (dims == 3)
        dz = connectivity(i,2);
      else
        dz = 0;

      lookup[ std::tuple<int, int, int>(dx,dy,dz) ] = i;
    }
  };

  double operator () (double x1,double y1,double z1, double x2, double y2, double z2)
  {
    int dx,dy,dz;

    dx = (int) x2 - x1;
    dy = (int) y2 - y1;

    if (dims == 3)
    {
      dz = (int) z2 - z1;
      int index = lookup[std::tuple<int, int, int>(dx,dy,dz)];
      return data(x1,y1,z1, index);
    }
    else
    {
      dz = 0;
      int index = lookup[std::tuple<int, int, int>(dx,dy,dz)];
      return data(x1,y1, index);
    }
  }

protected:
  std::map<std::tuple<int, int, int>, int> lookup;
  const matrix<double> data;
  const matrix<int> connectivity;
  unsigned char dims;
};