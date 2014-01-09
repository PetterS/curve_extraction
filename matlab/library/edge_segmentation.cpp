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

// Gives head and tail of the edge
Mesh::Point  tail_of_edge(int edge_num, const matrix<int>& connectivity)
{
  int num_points_per_element =  connectivity.M;
  int tail    = edge_num / num_points_per_element;

  return make_point(tail);
}

// Gives head and tail of the edge
Mesh::Point  head_of_edge(int edge_num, const matrix<int>& connectivity)
{
  int num_points_per_element =  connectivity.M;
  int edgeid = edge_num % num_points_per_element;

  Mesh::Point tail_point = tail_of_edge(edge_num, connectivity);

  return Mesh::Point( tail_point.x + connectivity(edgeid,0),
                      tail_point.y + connectivity(edgeid,1),
                      tail_point.z + connectivity(edgeid,2));
}


std::vector<Mesh::Point>  edgepath_to_points(const std::vector<int>& path, const matrix<int>& connectivity)
{
  std::vector<Mesh::Point> point_vector;
  Mesh::Point first();

  // Start point
  if (path.size() > 0) {
     Mesh::Point tail = tail_of_edge(path[0], connectivity);
     point_vector.push_back(tail);
  }

  for (int i = 0; i < path.size(); i++) {
     Mesh::Point head = head_of_edge(path[i], connectivity);
     point_vector.push_back(head);
  }

  return point_vector;
}



void edge_segmentation( std::vector<Mesh::Point>& points,
                        double& run_time,
                        int& evaluations,
                        double& cost,
                        const matrix<unsigned char>& mesh_map,
                        PieceWiseConstant& data_term,
                        const matrix<int>& connectivity,
                        const InstanceSettings& settings,
                        const PointSets& start_sets,
                        const PointSets& end_sets,
                        const std::vector<double>& voxeldimensions,
                        const ShortestPathOptions& options,
                        matrix<double>& visit_time)
{
  // Create functor handling regularization costs
  length_cost_functor length_cost(voxeldimensions, settings.length_penalty);
  curvature_cost_functor curvature_cost(voxeldimensions, settings.curvature_penalty, settings.curvature_power);

  // Some notation for the edge graph
  // Elements corresponds to points in the original graph
  // Points corresponds to edges in the original graph
  // Edges correspond to edgepairs in the original graph
  int num_elements = mesh_map.numel();
  int num_points_per_element = connectivity.M;
  int num_edges_per_point = num_points_per_element*num_points_per_element;

  // Total
  int num_edges = num_points_per_element*num_elements;

  // Filling the cache
  // connectivity.M is the number of edges from each  each node
  std::vector<double> regularization_cache(num_edges_per_point);

  int x,y,z, x2,y2,z2,x3,y3,z3, element_number, element_number_2;

  for (int i = 0; i < connectivity.M; i++) {
    x2 = connectivity(i,0);
    y2 = connectivity(i,1);
    z2 = connectivity(i,2);

  for (int j = 0; j < connectivity.M; j++) {
    x3 = x2 + connectivity(j,0);
    y3 = y2 + connectivity(j,1);
    z3 = z2 + connectivity(j,2);

    int n = i*num_points_per_element +j;

     regularization_cache[n] = curvature_cost(0,0,0,x2,y2,z2,x3,y3,z3)
                            + length_cost(x2,y2,z2,x3,y3,z3);
  }
  }

  if (verbose)
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

  // Extra start sets.
  for (int i = 0; i < start_sets.size(); ++i) {
    const auto points = start_sets[i];
    for (int j = 0; j < points.size() - 1; ++j)
    {
      int element_number = sub2ind(points[j]);

      int dx = points[j+1].x - points[j].x;
      int dy = points[j+1].y - points[j].y;
      int dz = points[j+1].z - points[j].z;

      int k;
      for (k = 0; k < connectivity.M; k++)
      {
        if (  (std::abs(dx - connectivity(k,0)) < 1e-6) &&
              (std::abs(dy - connectivity(k,1)) < 1e-6) &&
              (std::abs(dz - connectivity(k,2)) < 1e-6) )
          break;
      }

      if (k == connectivity.M)
        mexErrMsgTxt("Unable to find edge \n");

      start_set.insert(element_number*num_points_per_element + k);
    }
  }

  // Extra end sets.
  for (int i = 0; i < end_sets.size(); ++i) {
    const auto points = end_sets[i];
    for (int j = 0; j < points.size() - 1; ++j) {
      int element_number = sub2ind(points[j]);

      int dx = points[j+1].x - points[j].x;
      int dy = points[j+1].y - points[j].y;
      int dz = points[j+1].z - points[j].z;

      int k;
      for (k = 0; k < connectivity.M; k++)
      {
        if (  (std::abs(dx - connectivity(k,0)) < 1e-6) &&
              (std::abs(dy - connectivity(k,1)) < 1e-6) &&
              (std::abs(dz - connectivity(k,2)) < 1e-6) )
          break;
      }

      if (k == connectivity.M)
        mexErrMsgTxt("Unable to find edge \n");

      end_set.insert(element_number*num_points_per_element + k);
    }
  }

  int e_super = num_edges;
  std::set<int> super_edge;
  super_edge.insert(e_super);

  if (verbose)
    mexPrintf("done.\n");

  auto get_neighbors =
    [ &evaluations, &data_term, &num_points_per_element, &regularization_cache,
      &e_super, &start_set, &connectivity, &length_cost]
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

        int x1,y1,z1,
            x2,y2,z2;

        tie(x1,y1,z1) = ind2sub(root); 

        x2 = x1 + connectivity(k,0);
        y2 = y1 + connectivity(k,1);
        z2 = z1 + connectivity(k,2);

        float cost = data_term.evaluate_line_integral<double>
                      (x1, y1, z1,
                       x2, y2, z2);

        cost += length_cost(x1,y1,z1,x2,y2,z2);              

        neighbors->push_back(Neighbor(*itr, cost));
      }
    } else
    {
      // "Grow" out an edge pair passing to points p2 and p3.
      int root, edge_id_1;
      tie(root, edge_id_1) = root_and_edge(e, connectivity);

      int x1,y1,z1, x2,y2,z2;

      tie(x1,y1,z1) = ind2sub(root);

      x2 = x1 + connectivity(edge_id_1,0);
      y2 = y1 + connectivity(edge_id_1,1);
      z2 = z1 + connectivity(edge_id_1,2);

      if (validind(x2,y2,z2))
      {
        neighbors->resize(num_points_per_element);

        #ifdef USE_OPENMP
        #pragma omp parallel for
        #endif

        for (int edge_id_2 = 0; edge_id_2 < num_points_per_element; ++edge_id_2)
        {
          int x3 = x2 + connectivity(edge_id_2,0);
          int y3 = y2 + connectivity(edge_id_2,1);
          int z3 = z2 + connectivity(edge_id_2,2);

           // Id of destination, (element_number*num_points_per_element + edge_id)
          int dest = sub2ind(x2,y2,z2)*num_points_per_element + edge_id_2;
          float cost;

          if (validind(x3,y3,z3))
          {
            // Adjacency id  
            cost = data_term.evaluate_line_integral<double>
                    (x2, y2, z2,
                     x3, y3, z3);

            // Lookup id
            int edge_type =  edge_id_1*num_points_per_element + edge_id_2;
            cost += regularization_cache[edge_type];
          }

          else {
            cost = std::numeric_limits<float>::infinity();
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
    int p;
    tie(p, std::ignore) = root_and_edge(e, connectivity);
    return heuristic_options.distance[p];
  };

  std::function<double(int)>* lower_bound_pointer = nullptr;

  if (settings.use_a_star) {
    // Call node_segmentation. It solves the same problem, but without
    // the curvature term.
    // We tell it to compute all distances, and we will get a
    // vector of the distance from any node to the end set.
    // (node_segmentation switches the start and end sets.)
    double heuristic_runtime = 0;
    int heuristic_evaluations = 0;
    double heuristic_cost = 0;

    node_segmentation(points,
                      heuristic_runtime,
                      heuristic_evaluations,
                      heuristic_cost,
                      mesh_map,
                      data_term,
                      connectivity,
                      settings,
                      start_sets,
                      end_sets,
                      voxeldimensions,
                      heuristic_options,
                      visit_time);

    lower_bound_pointer = &lower_bound;
  }

  if (verbose)
    mexPrintf("Computing shortest curvature ...");

  std::vector<int> path_edges;
  double start_time = ::get_wtime();
  evaluations = 0;

  cost = shortest_path(num_edges+1,
                       super_edge,
                       end_set,
                       get_neighbors,
                       &path_edges,
                       lower_bound_pointer,
                       options);


  double end_time = ::get_wtime();
  run_time = end_time - start_time;

  path_edges.erase(path_edges.begin());  // Remove super edge
  points = edgepath_to_points(path_edges, connectivity);

  if (verbose)
  {
    mexPrintf("done.\n");
    mexPrintf("Running time:  %g (s), ", run_time);
    mexPrintf("Evaluations: %d, ", evaluations);
    mexPrintf("Path length: %d, ", path_edges.size() );
    mexPrintf("Cost:    %g. \n", cost);
  }

  // Store visit time
  if (options.store_visited) 
  {
    // Initialize.
    for (int i = 0; i < visit_time.numel(); ++i) 
        visit_time(i) = -1; 
 
    // Go through each each edge stored in visit time
    // if it has been visited then it's != -1
    std::vector<Mesh::Point> point_vector(2, make_point(0));
    for (int i = 0; i < options.visit_time.size(); i++)
    {
      if (options.visit_time[i] == -1)
        continue;

      point_vector[0] = tail_of_edge(i, connectivity);
      point_vector[1] = head_of_edge(i, connectivity);

      for (Mesh::Point p : point_vector)
      {
        if (!validind(p))
          continue;

        double visit_value = visit_time(p.x, p.y, p.z);

        if ( (visit_value == -1) || 
            ( (visit_value >= 0) && (visit_value > options.visit_time[i]) ) )
        {
          visit_time(p.x, p.y, p.z) = options.visit_time[i];
        }
      }
    }
  }
 }
