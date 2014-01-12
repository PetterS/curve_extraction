#include "curve_segmentation.h"

// Main nodes (start and end set flipped to accommodate for A*)
void node_segmentation(std::vector<Mesh::Point>& points,
                      double& run_time,
                      int& evaluations,
                      double& cost,
                      const matrix<unsigned char>& mesh_map,
                      PieceWiseConstant& data_term,
                      const matrix<int>& connectivity,
                      const InstanceSettings& settings,
                      const PointSets& start_sets,
                      const PointSets& end_sets,
                      const std::vector<double>& voxel_dimensions,
                      const ShortestPathOptions& options,
                      matrix<double>& visit_time,
                      matrix<int>& shortest_path_tree)
{
  // Create functor handling regularization costs
  length_cost_functor length_cost(voxel_dimensions, settings.length_penalty);
  
  // Precalculate regularization cost for every item connectivity
  std::vector<double> regularization_cache(connectivity.M);

  for (int k = 0; k < connectivity.M; k++)
  {
    int x = connectivity(k,0);
    int y = connectivity(k,1);
    int z = connectivity(k,2);

    regularization_cache[k] = length_cost(0,0,0,x,y,z);
  }

  auto get_neighbors =
    [&evaluations, &data_term, &connectivity, 
      &regularization_cache, &voxel_dimensions]
    (int n, std::vector<Neighbor>* neighbors) -> void
  {
    evaluations++;

    int x1,y1,z1;
    int x2,y2,z2;

    tie(x1,y1,z1) = ind2sub(n);

    if (validind(x1,y1,z1))
    {
      for (int k = 0; k < connectivity.M; k++)
      {
        x2 = x1 - connectivity(k,0);
        y2 = y1 - connectivity(k,1);
        z2 = z1 - connectivity(k,2);

        if (validind(x2,y2,z2))
        {
          // Unary
          float cost = data_term.evaluate_line_integral<double>
                      (x1,y1,z1,
                       x2,y2,z2);

          // Length reg;
          cost += regularization_cache[k];
          int dest = sub2ind(x2, y2, z2);

          neighbors->push_back(Neighbor(dest, cost));
        }
      }
    }
  };

  if (verbose)
    mexPrintf("Creating start/end sets without mesh...");

  // start and end set
  std::set<int> start_set, end_set;

  for (int i = 0; i < mesh_map.numel(); i++)
  {
      if ( mesh_map(i) == 3 )
          start_set.insert(i);

      if ( mesh_map(i) == 2 )
          end_set.insert(i);
  }

  // Extra start sets.
  for (int i = 0; i < start_sets.size(); ++i) {
    const auto& points = start_sets[i];
    for (int j = 0; j < points.size(); ++j) {

      int p = sub2ind(points[j]);
      end_set.insert(p);
    }
  }

  // Extra end sets.
  for (int i = 0; i < end_sets.size(); ++i) {
    const auto& points = end_sets[i];
    for (int j = 0; j < points.size(); ++j) {
      int p = sub2ind(points[j]);
      start_set.insert(p);
    }
  }

  if (verbose)
    mexPrintf("Computing shortest distance ...");

  std::vector<int> path_nodes;
  double start_time = ::get_wtime();
  evaluations = 0;

  //
  // NOTE!
  //
  // Length-based shortest path switches the start and end sets.
  // The code below is not a typo! It is equivalent for the best
  // path, but not for the distance map to the end set.
  //
  cost = shortest_path( mesh_map.numel(),
                        start_set,
                        end_set,
                        get_neighbors,
                        &path_nodes,
                        0,
                        options);
  std::reverse(path_nodes.begin(), path_nodes.end());

  double end_time = ::get_wtime();
  run_time = end_time - start_time;

  // Convert from inds to points
  for (auto itr = path_nodes.begin();
       itr != path_nodes.end();
       itr++)
  {
    points.push_back( make_point(*itr) );
  }

  if (verbose)
  {
    mexPrintf("done. \n");
    mexPrintf("Running time:  %g (s), ", run_time);
    mexPrintf("Evaluations: %d, ", evaluations);
    mexPrintf("Path length: %d, ", path_nodes.size() );
    mexPrintf("Cost:    %g. \n", cost);
  }

  if (options.store_visited) 
  {
    ASSERT(visit_time.numel() == options.visit_time.size());
    
    for (int i = 0; i < visit_time.numel(); i++)
      visit_time(i) = options.visit_time[i];
  }

  if (options.store_parents)
  {
    ASSERT(shortest_path_tree.numel() == options.parents.size());

    for (int i = 0; i < shortest_path_tree.numel(); i++)
      shortest_path_tree(i) = options.parents[i];
  }
}
