#include "curve_segmentation.h"

// Main nodes (start and end set flipped to accommodate for A*)
template<typename Data_cost, typename Pair_cost>
void node_segmentation( const matrix<double>& data,
                        const matrix<unsigned char>& mesh_map,
                        const matrix<int>& connectivity,
                        InstanceSettings& settings,
                        ShortestPathOptions& options,
                        SegmentationOutput& output
                      )
{
  bool reverse_direction = settings.use_a_star;

  Data_cost data_cost(data, connectivity, settings);
  Pair_cost pair_cost(data, settings);
  Delta_point delta_point(connectivity, reverse_direction);

  // Check overflow
  if (max_index < data.numel())
    mexErrMsgTxt("Problem is too large, index will overflow.");

  bool cacheable = true;
  if ( (pair_cost.data_dependent) && (settings.penalty[0] > 0) )
      cacheable = false;

  std::vector<double> regularization_cache(delta_point.size());

   // Pre-calculate regularization cost for every connectivity
  if (cacheable) 
  {
    Point p1 = make_point(0);

    for (int k = 0; k < delta_point.size(); k++)
    {
      Point p2 = delta_point(p1,k);

      if (!reverse_direction)
        regularization_cache[k] = pair_cost(p1.xyz, p2.xyz);
      else
        regularization_cache[k] = pair_cost(p2.xyz, p1.xyz);
    }
  }

  int evaluations = 0;
  auto get_neighbors =
    [&evaluations, &data_cost, 
      &regularization_cache, &cacheable, 
      &pair_cost, &delta_point, &reverse_direction]
    (int n, std::vector<Neighbor>* neighbors) -> void
  {
    evaluations++;
    Point p1 = make_point(n);

    neighbors->resize(delta_point.size());

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int k = 0; k < delta_point.size(); k++)
    {
      Point p2 = delta_point(p1,k);
      int dest;
      double cost;

      if (valid_point(p2))
      {
        dest = point2ind(p2);

        if (!reverse_direction)
        {
          cost = data_cost(p1.xyz, p2.xyz);
        
          if (cacheable)
            cost += regularization_cache[k];
          else
            cost += pair_cost(p1.xyz,p2.xyz);
        }
        else
        {
          cost = data_cost(p2.xyz, p1.xyz);

          if (cacheable)
            cost += regularization_cache[k];
          else
            cost += pair_cost(p2.xyz,p1.xyz);

        } 
      } 
      else {
        cost = std::numeric_limits<double>::infinity();
        dest = 0;
      }

      (*neighbors)[k] = Neighbor(dest, cost);
    }
  };

  if (settings.verbose)
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

  if (settings.verbose)
    mexPrintf("Computing shortest distance ...");

  std::vector<int> path_nodes;
  double start_time = ::get_wtime();
  evaluations = 0;

  if (reverse_direction)
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

  if (settings.verbose)
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