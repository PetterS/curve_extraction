#pragma once

// Explicitly defined cost for each edge.
class Edge_data_cost 
{
public:
  Edge_data_cost(
    const matrix<double>& data,
    const matrix<int>& connectivity,
    const InstanceSettings& settings
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

  template<typename R>
  R operator()(const R* const point1, const R* const point2) const
  {
    int dx,dy,dz;

    dx = (int) point2[0] - point1[0];
    dy = (int) point2[1] - point1[1];

    if (dims == 3)
    {
      dz = (int) point2[2] - point1[2];
      int index = lookup[std::tuple<int, int, int>(dx,dy,dz)];
      return data(point1[0], point1[1], point1[2], index);
    }
    else
    {
      dz = 0;
      int index = lookup[std::tuple<int, int, int>(dx,dy,dz)];
      return data(point1[0], point1[1], index);
    }
  }

protected:
  mutable std::map<std::tuple<int, int, int>, int> lookup;
  const matrix<double> data;
  const matrix<int> connectivity;
  unsigned char dims;
};