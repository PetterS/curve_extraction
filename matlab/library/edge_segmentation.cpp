#include "curve_segmentation.h"

// The indexing:
// Assume we have M neighbors in connectivity.
// Then a edge (e) starting in node i
// have index i*M + e
// The edge numbering is defined in the connectivity matrix.

std::tuple<int,int> root_and_edge(int edge_num, const matrix<int>& connectivity)
{
  int divisor =  connectivity.M;

  int root_node  = edge_num / divisor;
  int edge_id    = edge_num % divisor;

  return std::make_tuple(root_node,edge_id);
}

// Gives tail of the edge
Point  tail_of_edge(int edge_num, const matrix<int>& connectivity)
{
  int num_points_per_element =  connectivity.M;
  int tail    = edge_num / num_points_per_element;

  return make_point(tail);
}

// Gives head and tail of the edge
Point  head_of_edge(int edge_num, const matrix<int>& connectivity)
{
  int num_points_per_element =  connectivity.M;
  int edgeid = edge_num % num_points_per_element;

  Point tail_point = tail_of_edge(edge_num, connectivity);

  return Point( tail_point[0] + connectivity(edgeid,0),
                tail_point[1] + connectivity(edgeid,1),
                tail_point[2] + connectivity(edgeid,2));
}


std::vector<Point>  edgepath_to_points(const std::vector<int>& path, const matrix<int>& connectivity)
{
  std::vector<Point> point_vector;
  Point first();

  // Start point
  if (path.size() > 0) {
     Point tail = tail_of_edge(path[0], connectivity);
     point_vector.push_back(tail);
  }

  for (int i = 0; i < path.size(); i++) {
     Point head = head_of_edge(path[i], connectivity);
     point_vector.push_back(head);
  }

  return point_vector;
}

template<typename nodeT, typename edgeT>
void store_results_edge(matrix<nodeT>& node_container, std::vector<edgeT>& edge_container, const matrix<int>& connectivity)
{
  // Initialize.
  for (int i = 0; i < node_container.numel(); ++i)
      node_container(i) = -1;

  // Go through each each edge stored in visit time
  // if it has been visited then it's != -1
  std::vector<Point> point_vector(2, make_point(0));
  for (int i = 0; i < edge_container.size(); i++)
  {
    if (edge_container[i] == -1)
      continue;

    point_vector[0] = tail_of_edge(i, connectivity);
    point_vector[1] = head_of_edge(i, connectivity);

    for (Point p : point_vector)
    {
      if (!valid_point(p))
        continue;

      nodeT visit_value = node_container(p[0], p[1], p[2]);

      if ( (visit_value == -1) ||
          ( (visit_value >= 0) && (visit_value > edge_container[i]) ) )
      {
        node_container(p[0], p[1], p[2]) = edge_container[i];
      }
    }
  }

  return;
}



template<typename Data_cost, typename Length_cost, typename Curvature_cost>
void edge_segmentation( const matrix<double>& data,
                        const matrix<unsigned char>& mesh_map,
                        const matrix<int>& connectivity,
                        InstanceSettings& settings,
                        ShortestPathOptions& options,
                        SegmentationOutput& output)
{
  Data_cost data_cost(data, connectivity, settings.voxel_dimensions);
  Length_cost length_cost(data, settings.voxel_dimensions, settings.length_penalty);
  Curvature_cost curvature_cost(data, settings.voxel_dimensions, settings.curvature_penalty, settings.curvature_power);
  Delta_point delta_point(connectivity);

  bool cacheable = true;
  if ( (length_cost.data_depdent) && (settings.length_penalty > 0) )
      cacheable = false;

  if ( (curvature_cost.data_depdent) && (settings.curvature_penalty > 0) )
      cacheable = false;

  // Some notation for the edge graph
  // Elements corresponds to points in the original graph
  // Points corresponds to edges in the original graph
  // Edges correspond to edgepairs in the original graph
  int num_elements = mesh_map.numel();
  int num_points_per_element = connectivity.M;
  int num_edges_per_point = num_points_per_element*num_points_per_element;

  if (max_index < mesh_map.numel()*num_points_per_element)
      mexErrMsgTxt("Problem is too large, index will overflow. Try to remove curvature penalty.");

  // Total
  int num_edges = num_points_per_element*num_elements;

  // Filling the cache
  // connectivity.M is the number of edges from each  each node
  std::vector<double> regularization_cache(num_edges_per_point);

  int x,y,z, x2,y2,z2, element_number, element_number_2;
  if (cacheable)
  {
    Point p1 = make_point(0);

    for (int i = 0; i < connectivity.M; i++) {
      Point p2 = delta_point(p1,i);

    for (int j = 0; j < connectivity.M; j++) {
      Point p3 = delta_point(p2,j);
      int n = i*num_points_per_element +j;

      regularization_cache[n] = curvature_cost( p1.xyz, p2.xyz, p3.xyz)
                              + length_cost(            p2.xyz, p3.xyz);
    }
    }
  }

  if (settings.verbose)
    mexPrintf("Creating start/end sets...");

  // start and end set
  std::set<int> start_set, end_set;

  // Add edges according to mesh_map
  for (x = 0; x < M; x++) {
  for (y = 0; y < N; y++) {
  for (z = 0; z < O; z++) {

    // Add all edges
    if (mesh_map(x,y,z) ==2)
    {
        element_number = sub2ind(x,y,z);

        for (int k = 0; k < num_points_per_element; k++)
        {
          if ( validind(x + connectivity(k,0),
                        y + connectivity(k,1),
                        z + connectivity(k,2)))
          {
              start_set.insert(element_number*num_points_per_element + k);
          }
        }
    }

    // End set, check which edges goes into this one
    if (mesh_map(x,y,z) == 3)
    {
        for (int k = 0; k < connectivity.M; k++)
        {
          x2 = x - connectivity(k,0);
          y2 = y - connectivity(k,1);
          z2 = z - connectivity(k,2);

          // The edge with head at n has root at n2
          // it's edge number k.
          if ( validind(x2,y2,z2) )
          {
            element_number_2 = sub2ind(x2,y2,z2);
            end_set.insert(element_number_2*num_points_per_element + k);
          }
        }
      }
  }
  }
  }

  int e_super = num_edges;
  std::set<int> super_edge;
  super_edge.insert(e_super);

  if (settings.verbose)
    mexPrintf("done.\n");

  int evaluations = 0;
  auto get_neighbors =
    [ &evaluations, &data_cost, &num_points_per_element, &regularization_cache,
      &e_super, &start_set, &connectivity, &length_cost, 
      &cacheable, &curvature_cost, &delta_point]
    (int e, std::vector<Neighbor>* neighbors) -> void
  {
    evaluations++;

    if (e == e_super)
    {
      for (auto itr = start_set.begin();
           itr != start_set.end();
           itr++)
      {
        int root, k;
        tie(root, k) = root_and_edge(*itr, connectivity);

        Point p1 = make_point(root);
        Point p2 = delta_point(p1, k);

        double cost  = data_cost(   p1.xyz, p2.xyz);
        cost        += length_cost( p1.xyz, p2.xyz);

        neighbors->push_back(Neighbor(*itr, cost));
      }
    } else
    {
      // "Grow" out an edge pair passing to points p2 and p3.
      int root, edge_id_1;
      tie(root, edge_id_1) = root_and_edge(e, connectivity);

      Point p1 = make_point(root);
      Point p2 = delta_point(p1, edge_id_1);
      
      if (valid_point(p2))
      {
        neighbors->resize(num_points_per_element);

        #ifdef USE_OPENMP
        #pragma omp parallel for
        #endif

        for (int edge_id_2 = 0; edge_id_2 < num_points_per_element; ++edge_id_2)
        {
          Point p3 = delta_point(p2, edge_id_2);

           // Id of destination, (element_number*num_points_per_element + edge_id)
          int dest = point2ind(p2)*num_points_per_element + edge_id_2;
          double cost;

          if (valid_point(p3))
          {
            // Adjacency id
            cost = data_cost(p2.xyz, p3.xyz);

            if (cacheable)
            {
              // Lookup id
              int edge_type =  edge_id_1*num_points_per_element + edge_id_2;
              cost += regularization_cache[edge_type];
            } else
            {
              cost += curvature_cost(p1.xyz,p2.xyz, p3.xyz);
              cost += length_cost   (       p2.xyz, p3.xyz);
            }
          }

          else {
            cost = std::numeric_limits<double>::infinity();
          }

          (*neighbors)[edge_id_2] = Neighbor(dest, cost);
        }
      }
    }
  };

  ShortestPathOptions heuristic_options;
  heuristic_options.compute_all_distances = true;

  // The lower bound function is just the distance
  // without curvature taken into account.
  std::function<double(int)> lower_bound =
    [&heuristic_options, &connectivity]
    (int e) -> double
  {
    Point p = head_of_edge(e, connectivity);
    return heuristic_options.distance[point2ind(p)];
  };

  std::function<double(int)>* lower_bound_pointer = nullptr;

  if (settings.use_a_star && !options.store_parents) {
    // Call node_segmentation. It solves the same problem, but without
    // the curvature term.
    // We tell it to compute all distances, and we will get a
    // vector of the distance from any node to the end set.
    // (node_segmentation switches the start and end sets.)
    // If all parents are to be stored A* will not help.
    double heuristic_runtime = 0;
    int heuristic_evaluations = 0;
    double heuristic_cost = 0;
    matrix<int> empty_matrix;
    matrix<double> empty_double_matrix;

    SegmentationOutput heuristic_output
    (output.points, heuristic_runtime, heuristic_evaluations, heuristic_cost, output.visit_time, empty_matrix, empty_double_matrix);

    node_segmentation<Data_cost, Length_cost>
                     (data,
                      mesh_map,
                      connectivity,
                      settings,
                      heuristic_options,
                      heuristic_output);

    lower_bound_pointer = &lower_bound;
  }

  if (settings.verbose)
    mexPrintf("Computing shortest curvature ...");

  std::vector<int> path_edges;
  double start_time = ::get_wtime();
  evaluations = 0;

  // This can be a bit confusing: when we use line-graphs
  // the parent is built into the edge so we do need to store them
  // we DO however need to store the visit order in order to resolve conflicts.
  if (options.store_parents)
    options.store_visited = true;

  options.store_parents = false;

  output.cost = shortest_path(num_edges+1,
                       super_edge,
                       end_set,
                       get_neighbors,
                       &path_edges,
                       lower_bound_pointer,
                       options);

  // Code clarity
  options.store_parents = settings.store_parents;

  double end_time = ::get_wtime();
  output.run_time = end_time - start_time;

  path_edges.erase(path_edges.begin());  // Remove super edge
  output.points = edgepath_to_points(path_edges, connectivity);

  output.evaluations = evaluations;
  if (settings.verbose)
  {
    mexPrintf("done.\n");
    mexPrintf("Running time:  %g (s), ", output.run_time);
    mexPrintf("Evaluations: %d, ", output.evaluations);
    mexPrintf("Path length: %d, ", path_edges.size() );
    mexPrintf("Cost:    %g. \n", output.cost);
  }

  if (settings.store_distances)
    store_results_edge<double,float>(output.distances, options.distance, connectivity);

  // Store visit time
  if (options.store_visited)
    store_results_edge<int,int>(output.visit_time, options.visit_time, connectivity);
 
  // Store parents
  // Conflicts are resolved by first visit.
  if (options.store_parents)
  {
    ASSERT(options.store_visited);

    // Initialize.
    for (int i = 0; i < output.shortest_path_tree.numel(); ++i)
        output.shortest_path_tree(i) = -1;

    // Go through each edge stored in visit time
    std::vector<Point> point_vector(2, make_point(0));
    for (int i = 0; i < options.visit_time.size(); i++)
    {
      point_vector[0] = tail_of_edge(i, connectivity);
      if (!valid_point(point_vector[0]))
        continue;

      point_vector[1] = head_of_edge(i, connectivity);
      if (!valid_point(point_vector[1]))
        continue;

      int time = output.visit_time( point_vector[1][0],
                                    point_vector[1][1],
                                    point_vector[1][2]);

      // Is this the edge which was here first?
      if (time == options.visit_time[i])
      {
        output.shortest_path_tree(  point_vector[1][0],
                                    point_vector[1][1],
                                    point_vector[1][2]) = point2ind(point_vector[0]);
      }
    }
  }
}