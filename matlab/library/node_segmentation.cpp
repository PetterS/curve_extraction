#include "curve_segmentation.h"

// Main nodes (start and end set flipped to accommodate for A*)
template<typename Data_cost, typename Length_cost>
void node_segmentation( const matrix<double>& data,
                        const matrix<unsigned char>& mesh_map,
                        const matrix<int>& connectivity,
                        InstanceSettings& settings,
                        ShortestPathOptions& options,
                        SegmentationOutput& output
                      )
{
  Data_cost data_cost(data, connectivity, settings.voxel_dimensions);
  Length_cost length_cost(data, settings.voxel_dimensions, settings.length_penalty);
  Delta_point delta_point(connectivity);

  bool cacheable = true;
  if ( (length_cost.data_depdent) && (settings.length_penalty > 0) )
      cacheable = false;

  std::vector<double> regularization_cache(connectivity.M);

   // Pre-calculate regularization cost for every connectivity
  if (cacheable) 
  {
    Point p1 = make_point(0);

    for (int k = 0; k < connectivity.M; k++)
    {
      Point p2 = delta_point(p1,k);
      regularization_cache[k] = length_cost(p1.xyz, p2.xyz);
    }
  }

  int evaluations = 0;
  auto get_neighbors =
    [&evaluations, &data_cost, &connectivity, 
      &regularization_cache, &cacheable, 
      &length_cost, &delta_point]
    (int n, std::vector<Neighbor>* neighbors) -> void
  {
    evaluations++;
    Point p1 = make_point(n);

    if (valid_point(p1))
    {
      for (int k = 0; k < connectivity.M; k++)
      {
        Point p2 = delta_point(p1,k);
        if (valid_point(p2))
        {
          // Unary
          double cost = data_cost(p1.xyz, p2.xyz);

          // Length reg;
          if (cacheable)
            cost += regularization_cache[k];
          else
            cost += length_cost(p1.xyz,p2.xyz);

          int dest = point2ind(p2);

          neighbors->push_back(Neighbor(dest, cost));
        }
      }
    }
  };

  if (verbose)
    mexPrintf("Creating start/end sets without mesh...");

  // start and end set
  std::set<int> end_set, start_set;

  for (int i = 0; i < mesh_map.numel(); i++)
  {
      if ( mesh_map(i) == 3 )
          end_set.insert(i);

      if ( mesh_map(i) == 2 )
          start_set.insert(i);
  }

  if (verbose)
    mexPrintf("Computing shortest distance ...");

  std::vector<int> path_nodes;
  double start_time = ::get_wtime();
  evaluations = 0;

  if (settings.use_a_star)
  {
    //
    // NOTE!
    //
    // Switching the start and end sets.
    // The code below is not a typo! It is equivalent for the best
    // path, but not for the distance map to the end set.
    //
    output.cost = shortest_path( mesh_map.numel(),
                          end_set,
                          start_set, 
                          get_neighbors,
                          &path_nodes,
                          0,
                          options);

    std::reverse(path_nodes.begin(), path_nodes.end());
  } else
  {
    output.cost = shortest_path( mesh_map.numel(),
                          start_set,
                          end_set, 
                          get_neighbors,
                          &path_nodes,
                          0,
                          options);
  }

  double end_time = ::get_wtime();
  output.run_time = end_time - start_time;

  // Convert from inds to points
  for (auto id : path_nodes)
    output.points.push_back( make_point(id) );

  output.evaluations = evaluations;

  if (verbose)
  {
    mexPrintf("done. \n");
    mexPrintf("Running time:  %g (s), ", output.run_time);
    mexPrintf("Evaluations: %d, ", output.evaluations);
    mexPrintf("Path length: %d, ", path_nodes.size() );
    mexPrintf("Cost:    %g. \n", output.cost);
  }

  if (settings.store_distances)
  {
    for (int i = 0; i < output.distances.numel(); i++)
      output.distances(i) = options.distance[i];
  }

  if (options.store_visited) 
  {
    for (int i = 0; i < output.visit_time.numel(); i++)
      output.visit_time(i) = options.visit_time[i];
  }

  if (options.store_parents)
  {
    for (int i = 0; i < output.shortest_path_tree.numel(); i++)
      output.shortest_path_tree(i) = options.parents[i];
  }  
}