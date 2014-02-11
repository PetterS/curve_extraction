#include "curve_segmentation.h"
#define EDGEPAIR_SEGMENTATION

// The indexing:
// Assume we have M neighbors in connectivity.
// Then a pair edges starting in node i with edge e1 and e2
// have index i*M*M + e*M + e2;
// where e1,e2 are indices defined by the connectivity matrix.

std::tuple<int,int> decompose_edgepair(int edge_num, const matrix<int>& connectivity)
{
  int divisor =  connectivity.M;

  int root_node  = edge_num / divisor;
  int edge_id    = edge_num % divisor;

  return std::make_tuple(root_node,edge_id);
}

std::tuple<int,int> 
decompose_pair_of_edgepairs(int edgepair_num, const matrix<int>& connectivity)
{
  int divisor = connectivity.M*connectivity.M;

  int root_node   = edgepair_num / divisor;
  int edgepair_id = edgepair_num % divisor;

  return std::make_tuple(root_node, edgepair_id);
}

// Given edgeapir_id return the three element id's associated with that edgepair 
std::tuple<int, int, int> 
points_in_a_edgepair(int edgepair_num, const matrix<int>& connectivity)
{
  int root, edgepair_id;
  tie(root, edgepair_id) = decompose_pair_of_edgepairs(edgepair_num, connectivity);

  int e1, e2;
  tie(e1, e2) = decompose_edgepair(edgepair_id, connectivity);

  int x,y,z,x2,y2,z2, x3,y3,z3;
  tie(x,y,z)  = ind2sub(root);

  x2 = x + connectivity(e1,0);
  y2 = y + connectivity(e1,1);
  z2 = z + connectivity(e1,2); 

  x3 = x2 + connectivity(e2,0);
  y3 = y2 + connectivity(e2,1);
  z3 = z2 + connectivity(e2,2); 

  int q2,q3;
  q2 = sub2ind(x2,y2,z2);
  q3 = sub2ind(x3,y3,z3);

  return std::make_tuple(root,q2,q3);
}


std::vector<Mesh::Point> pairpath_to_points(const std::vector<int>& path, const matrix<int>& connectivity)
{
  std::vector<Mesh::Point> point_vector;

  // Start points
  if (path.size() > 0) {
    int p1, p2;
    tie(p1, p2, ignore) = points_in_a_edgepair(path[0], connectivity);
    point_vector.push_back( make_point(p1) );
    point_vector.push_back( make_point(p2) );
  }

  for (int i = 0; i < path.size(); i++)
  {
    int p3;
    tie(ignore,ignore, p3) = points_in_a_edgepair(path[i], connectivity);
    point_vector.push_back( make_point(p3) );
  }
  
  return point_vector;
}

template<typename Data_cost, typename Length_cost, typename Curvature_cost, typename Torsion_cost>
void  edgepair_segmentation(  const matrix<double>& data,
                              const matrix<unsigned char>& mesh_map,
                              const matrix<int>& connectivity,
                              InstanceSettings& settings,
                              ShortestPathOptions& options,
                              SegmentationOutput& output
                             )
{  
  Data_cost data_cost(data, connectivity, settings.voxel_dimensions);
  Length_cost length_cost(data,settings.voxel_dimensions, settings.length_penalty);
  Curvature_cost curvature_cost(data, settings.voxel_dimensions, settings.curvature_penalty, settings.curvature_power);
  Torsion_cost torsion_cost(data, settings.voxel_dimensions, settings.torsion_penalty, settings.torsion_power);

  // Some notation for the edge graph
  // Elements corresponds to points in the original graph
  // Points corresponds to edge pairs  in the original graph
  // Edges correspond to pair of edgepairs in the original graph

  int num_elements = mesh_map.numel();
  int num_points_per_element = connectivity.M*connectivity.M;

  // Total
  int num_edges = num_points_per_element*num_elements;

  // Read mesh_map to find end and start set.
  std::set<int> start_set_pairs, end_set_pairs;

  int x1,y1,z1,
      x2,y2,z2,
      x3,y3,z3;

  for (int n = 0; n < mesh_map.numel(); ++n)
  {
    tie(x1,y1,z1) = ind2sub(n);

    for (int e1 = 0; e1 < connectivity.M; e1++) 
    { 
      x2 = x1 + connectivity(e1,0);
      y2 = y1 + connectivity(e1,1);
      z2 = z1 + connectivity(e1,2);

      if (!validind(x2,y2,z2))
        continue;

      for (int e2 = 0; e2 < connectivity.M; e2++) 
      {
        x3 = x2 + connectivity(e2,0);
        y3 = y2 + connectivity(e2,1);
        z3 = z2 + connectivity(e2,2);

        // Symmetric neighborhood leads to useless pair going back to itself.
        if ((x1 == x3) & (y1 == y3) & (z1 == z3))
          continue;

        if (!validind(x3,y3,z3))
           continue;

        int pair_id = sub2ind(x1,y1,z1)*(connectivity.M*connectivity.M) 
                 + connectivity.M*e1 + e2;

        if ( mesh_map( x1, y1, z1) == 2)
          start_set_pairs.insert(pair_id);

        if ( mesh_map( x3, y3, z3) == 3)
          end_set_pairs.insert(pair_id);
      }
    }
  }


  // Super_edge which all edges goes out from. This is needed
  // because otherwise the first edge will have 0 regularization
  int e_super = num_edges;
  std::set<int> super_edge;
  super_edge.insert(e_super);

  int evaluations = 0;
 auto get_neighbors_torsion =
    [&evaluations, &data_cost,
     &e_super, &connectivity, &start_set_pairs,
     &length_cost, &curvature_cost, &torsion_cost]
    (int ep, std::vector<Neighbor>* neighbors) -> void
  {
    evaluations++;

    // Special case
    if (ep == e_super)
    {
      for (auto itr = start_set_pairs.begin();
           itr != start_set_pairs.end();
           itr++)
      {
        int q1, q2, q3;
        tie(q1, q2, q3) = points_in_a_edgepair(*itr, connectivity);

        int x2,y2,z2,
            x3,y3,z3, 
            x4,y4,z4;

        tie(x2,y2,z2) = ind2sub(q1);
        tie(x3,y3,z3) = ind2sub(q2);
        tie(x4,y4,z4) = ind2sub(q3);

        double cost = data_cost(x2, y2, z2, x3, y3, z3);
              cost += data_cost(x3, y3, z3, x4, y4, z4);

        cost += curvature_cost(x2,y2,z2, x3,y3,z3, x4,y4,z4);
        cost += length_cost(x2,y2,z2,x3,y3,z3);
        cost += length_cost(x3,y3,z3,x4,y4,z4);

        neighbors->push_back(Neighbor(*itr, cost));
      }
    }

    // The usual case
    else
    {
      // Given an edgepair find adjacent pairs
      //
      // From "this nodes" all neighboring edgepairs start.
      // o -- o -- o -- o (edge pair)
      //      ^ this node.
      int root, edgepair_id;
      tie(root, edgepair_id) = decompose_pair_of_edgepairs(ep, connectivity);

      int e1, e2;
      tie(e1,e2) = decompose_edgepair(edgepair_id, connectivity);

      int x1,y1,z1, 
          x2,y2,z2,
          x3,y3,z3,
          x4,y4,z4;

      tie(x1,y1,z1) = ind2sub(root);

      x2 = x1 + connectivity(e1,0);
      y2 = y1 + connectivity(e1,1);
      z2 = z1 + connectivity(e1,2);

      x3 = x2 + connectivity(e2,0);
      y3 = y2 + connectivity(e2,1),
      z3 = z2 + connectivity(e2,2);

      for (int e3 = 0; e3 < connectivity.M; e3++) 
      {
        x4 = x3 + connectivity(e3,0);
        y4 = y3 + connectivity(e3,1),
        z4 = z3 + connectivity(e3,2);

        double cost;
        if (!validind(x4,y4,z4))
          continue;

        if ((x2 == x4) & (y2 == y4) & (z2 == z4))
          continue;

        // Unary cost
        cost = data_cost(x3,y3,z3, x4,y4,z4);

        cost += length_cost(x3,y3,z3,x4,y4,z4);
        cost += torsion_cost(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4);
        cost += curvature_cost(x2,y2,z2, x3,y3,z3, x4,y4,z4);

        // Destination id
        int edge_pair_id = connectivity.M*e2 + e3;

        // Index of neighboring edgepair.
        int dest = sub2ind(x2,y2,z2)*(connectivity.M*connectivity.M) + edge_pair_id;
        neighbors->push_back(Neighbor(dest, cost));
      }
    }
  };

  // Compute shortest path
  evaluations = 0;
  std::vector<int> path_pairs;

  if (options.store_parents)
    options.store_visited = true;
  options.store_parents = false;

  if (verbose)
    mexPrintf("Computing shortest distance ...");

  double start_time = ::get_wtime();
  output.cost = shortest_path( num_edges+1,
                        super_edge,
                        end_set_pairs,
                        get_neighbors_torsion,
                        &path_pairs,
                        0,
                        options);

  // Code clarity
  options.store_parents = settings.store_parents;

  double end_time = ::get_wtime();
  output.run_time = end_time - start_time;
  path_pairs.erase(path_pairs.begin()); // Remove super edge
  output.points = pairpath_to_points(path_pairs, connectivity);

  output.evaluations = evaluations;
  if (verbose)
  {
    mexPrintf("done. \n");
    mexPrintf("Running time:  %g seconds,", output.run_time);
    mexPrintf("Evaluations: %d,", output.evaluations);
    mexPrintf("Path length: %d,", path_pairs.size() );
    mexPrintf("Cost:    %g. \n", output.cost);
  }

  // Store visit time
  if (options.store_visited) 
  {
    // Initialize.
    for (int i = 0; i < output.visit_time.numel(); ++i) 
        output.visit_time(i) = -1; 
 
    // Go through each each edge stored in visit time
    // if it has been visited then it's != -1
    int p0,p1,p2;

    // No empty constructor.
    std::vector<Mesh::Point> point_vector(3, make_point(0));

    for (int i = 0; i < options.visit_time.size(); i++)
    {
      if (options.visit_time[i] == -1)
        continue;

      tie(p0,p1,p2) = points_in_a_edgepair(i, connectivity);
      
      point_vector[0] = make_point(p0);
      point_vector[1] = make_point(p1);
      point_vector[2] = make_point(p2);

      for (Mesh::Point p : point_vector)
      {
        if (!validind(p))
          continue;

        double visit_value = output.visit_time(p.x, p.y, p.z);

        if ( (visit_value == -1) || 
            ( (visit_value >= 0) && (visit_value > options.visit_time[i]) ) )
        {
          output.visit_time(p.x, p.y, p.z) = options.visit_time[i];
        }
      }
    }
  }  

    // Store parents
  // Conflicts are resolved by first visit.
  if (options.store_parents)
  {
    ASSERT(options.store_visited);

    // Initialize.
    for (int i = 0; i < output.shortest_path_tree.numel(); ++i)
        output.shortest_path_tree(i) = -1;

    // Go through each each edge stored in visit time
    // if it has been visited then it's != -1
    int p0,p1,p2;
    std::vector<Mesh::Point> point_vector(3, make_point(0));
    for (int i = 0; i < options.visit_time.size(); i++)
    {
      tie(p0,p1,p2) = points_in_a_edgepair(i, connectivity);
      point_vector[0] = make_point(p0);
      point_vector[1] = make_point(p1);
      point_vector[2] = make_point(p2);

      if (!validind(point_vector[0]))
        continue;

      if (!validind(point_vector[1]))
        continue;

      if (!validind(point_vector[2]))
        continue;

      // First edge
      int time = output.visit_time(point_vector[1].x,
                            point_vector[1].y,
                            point_vector[1].z);

      // Is this the edge which was here first?
      if (time == options.visit_time[i])
      {
        output.shortest_path_tree(point_vector[1].x,
                           point_vector[1].y,
                           point_vector[1].z) = sub2ind(point_vector[0]);
      }

      // Second edge
      time = output.visit_time(point_vector[2].x,
                        point_vector[2].y,
                        point_vector[2].z);

      // Is this the edge which was here first?
      if (time == options.visit_time[i])
      {
        output.shortest_path_tree(point_vector[2].x,
                           point_vector[2].y,
                           point_vector[2].z) = sub2ind(point_vector[1]);
      }
    }
  }
}