// Petter Strandmark 2013.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <tuple>
using std::ignore;
using std::tie;

#ifdef USE_OPENMP
#include <omp.h>
double get_wtime()
{
	return ::omp_get_wtime();
}
#else
#include <ctime>
double get_wtime()
{
	return std::time(0);
}
#endif

#include <curve_extraction/curvature.h>
#include <curve_extraction/mesh.h>
#include <curve_extraction/shortest_path.h>

using namespace curve_extraction;

void compare_dijkstras_and_astar(
	int n,
	const std::function<void(int, std::vector<Neighbor>*)>& neighbors,
	const std::function<double(int)>& heuristic,
	const std::set<int>& start_set,
	const std::set<int>& end_set,
	std::vector<int>* path,
	int* evaluations,
	int* heuristic_evaluations)
{
	ShortestPathOptions options;
	options.print_progress = true;

	// Compute the shortest path.
	std::cerr << "Computing Dijkstras... ";
	double start_time = ::get_wtime();
	*evaluations = 0;
	double dist = shortest_path(n, start_set, end_set,
	                            neighbors, path, 0, options);
	double end_time = ::get_wtime();
	std::cerr << "done.\n";
	std::cerr << "Dijkstra's time: " << end_time - start_time << " seconds ";
	std::cerr << "(" << *evaluations << " evaluations).\n";

	std::cerr << "Computing A*... ";
	start_time = ::get_wtime();
	*evaluations = 0;
	*heuristic_evaluations = 0;
	double dist2 = shortest_path(n, start_set, end_set,
	                             neighbors, path, heuristic,
	                             options);
	end_time = ::get_wtime();
	std::cerr << "done.\n";
	std::cerr << "A* time        : " << end_time - start_time << " seconds ";
	std::cerr << "(" << *evaluations << " evaluations, ";
	std::cerr << *heuristic_evaluations << " heuristics).\n";

	if (std::abs(dist - dist2) > 1e-10) {
		std::cerr << dist << "  " << dist2 << std::endl;
		throw std::runtime_error("A* does not match Dijkstra's.");
	}

	std::cerr << "Shortest distance: " << dist << std::endl;
	std::cerr << "Path length: " << path->size() << std::endl;
	std::cerr << std::endl;
}

int main_function()
{
	const int n = 40;
	bool use_pairs = true;

	// Create the mesh and add all n*n points.
	Mesh mesh;
	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < n; ++y) {
			mesh.add_point(float(x), float(y), 0.0f);
		}
	}

	// This section defines the edges that are *not* part of the mesh.
	// All edges crossing a set of horizontal line segments are excluded.
	//
	// The line segments are  { y = y0
	//                        { x_min <= x <= x_max,
	// for y0 in {10.9, 10.1, 20.9, 20.1, 30.9, 30.1},
	//  x_min in {   0,    0, 2n/3,  2n/3,   0,    0} and
	//  x_max in { n/3,  n/3,    n,     n, n/3,  n/3}.
	std::vector<float> yvals;
	std::vector<float> xmin;
	std::vector<float> xmax;

	yvals.push_back(10.9f);
	xmin.push_back(-1);
	xmax.push_back(n / 3 - 0.1f);
	yvals.push_back(10.1f);
	xmin.push_back(-1);
	xmax.push_back(n / 3 - 0.1f);


	yvals.push_back(20.9f);
	xmin.push_back(2 * n / 3 + 0.1f);
	xmax.push_back(10000);
	yvals.push_back(20.1f);
	xmin.push_back(2 * n / 3 + 0.1f);
	xmax.push_back(10000);

	yvals.push_back(30.9f);
	xmin.push_back(-1);
	xmax.push_back(n / 3 - 0.1f);
	yvals.push_back(30.1f);
	xmin.push_back(-1);
	xmax.push_back(n / 3 - 0.1f);

	// This function returns whether a given edge should be
	// in the mesh.
	auto forbidden_edge = [n, &yvals, &xmin, &xmax]
	                      (float x1, float y1, float z1,
	                       float x2, float y2, float z2)
	                      -> bool
	{
		for (int i = 0; i < yvals.size(); ++i) {
			float y0 = yvals[i];

			float dy = y2 - y1;
			if (std::abs(dy) > 0) {
				// y1 + t * dy = y0  <=>
				// t = (y0 - y1) / dy
				float t = (y0 - y1) / dy;
				if (0.0 < t && t < 1.0) {
					float dx  = x2 - x1;
					float x0 = x1 + t * dx;
					if (x0 > xmax[i] && std::abs(dx) > 0.05) {
						return true;
					}
					if (x0 > xmax[i] + 0.2) {
						return true;
					}
					if (x0 < xmin[i] && std::abs(dx) > 0.05) {
						return true;
					}
					if (x0 < xmin[i] - 0.2) {
						return true;
					}
				}
			}
			else {
				if (std::abs(y1 - y0) < 0.05f &&
					x1 > xmax[i]) {
					return true;
				}
				if (std::abs(y2 - y0) < 0.05f &&
					x2 > xmax[i]) {
					return true;
				}
				if (std::abs(y1 - y0) < 0.05f &&
					x1 < xmin[i]) {
					return true;
				}
				if (std::abs(y2 - y0) < 0.05f &&
					x2 < xmin[i]) {
					return true;
				}
			}
		}
		return false;
	};

	// Add all allowed edges of length less than or
	// equal to four.
	double start_time = ::get_wtime();
	mesh.add_edges(4.0, forbidden_edge);
	double end_time = ::get_wtime();
	std::cerr << "Adding edges time: " << end_time - start_time << " seconds.\n";

	std::cerr << mesh.number_of_points() << " points and "
	          << mesh.number_of_edges() << " edges.\n";

	start_time = ::get_wtime();
	mesh.finish(use_pairs);
	end_time = ::get_wtime();
	std::cerr << "Mesh finalization time: " << end_time - start_time << " seconds." << std::endl;

	std::cerr << mesh.number_of_edge_pairs() << " edge pairs.\n\n";

	// Main function for computing the shortest path. Given an edge e, it
	// fills the neighbors vector with a list of the neighbors and their
	// distances.
	int evaluations = 0;
	auto get_neighbors =
		[&mesh, &evaluations]
		(int e, std::vector<Neighbor>* neighbors) -> void
	{
		evaluations++;

		static std::vector<int> adjacent;
		mesh.get_adjacent_edges(e, &adjacent);
		for (auto itr = adjacent.begin(); itr != adjacent.end(); ++itr) {
			float cost = mesh.edge_length(*itr);
			neighbors->push_back(Neighbor(*itr, cost));
		}
	};

	int heuristic_evaluations = 0;
	auto heuristic_lambda =
		[n, &mesh, &heuristic_evaluations]
		(int i) -> double
	{
		heuristic_evaluations++;

		int p1 = mesh.get_edge(i).first;
		int p2 = mesh.get_edge(i).second;
		float maxy = std::max(mesh.get_point(p1).y, mesh.get_point(p2).y);
		return std::max(0.0f, n - 1.5f - maxy);
	};
	std::function<double(int)> heuristic(heuristic_lambda);

	// The two sets between which to compute the shortest path.
	std::set<int> start_set, end_set;
	// Create the start and end set.
	for (int e = 0; e < mesh.number_of_edges(); ++e) {
		float y1 = mesh.get_point(mesh.get_edge(e).first).y;
		float y2 = mesh.get_point(mesh.get_edge(e).second).y;
		// Require that both y-coordinates are along the start line
		// because the first edge in the path is not counted.
		if (y1 < 0.5 && y2 < 0.5) {
			start_set.insert(e);
		}
		if (y1 > n - 1.5 || y2 > n - 1.5) {
			end_set.insert(e);
		}
	}

	// This vector will contain the shortest path.
	std::vector<int> path;

	// Perform the edge comparision
	std::cerr << "Edges\n";
	std::cerr << "-----------\n";
	compare_dijkstras_and_astar(mesh.number_of_edges(), get_neighbors,
	                            heuristic, start_set, end_set, &path,
								&evaluations, &heuristic_evaluations);

	// Function to specify the style of each edge in the SVG file
	// displaying the mesh.
	auto edge_style = [&mesh, &forbidden_edge](int e) -> std::string
	{
		return "stroke-width:0.01;stroke:#4f4f4f;";
	};

	mesh.start_SVG("2D.svg", edge_style);
	mesh.draw_path(path, "00aa00");

	// Neighborhood function for curvature computation.
	// The neighborhood structure is the same, but the costs
	// are different.
	auto get_neighbors_curvature =
		[&mesh, &evaluations]
		(int e, std::vector<Neighbor>* neighbors) -> void
	{
		evaluations++;

		static std::vector<int> adjacent;
		mesh.get_adjacent_edges(e, &adjacent);
		for (auto itr = adjacent.begin(); itr != adjacent.end(); ++itr) {
			int q1 = mesh.get_edge(e).first;
			int q2 = mesh.get_edge(e).second;
			int q3 = mesh.get_edge(*itr).first;
			int q4 = mesh.get_edge(*itr).second;
			int p1, p2, p3;
			if (q2 == q3) {
				p1 = q1;
				p2 = q2;
				p3 = q4;
			}
			else {
				throw std::runtime_error("Mesh error.");
			}

			float x1 = mesh.get_point(p1).x;
			float y1 = mesh.get_point(p1).y;
			float x2 = mesh.get_point(p2).x;
			float y2 = mesh.get_point(p2).y;
			float x3 = mesh.get_point(p3).x;
			float y3 = mesh.get_point(p3).y;

			float cost = compute_curvature<float>(x1,y1,0, x2,y2,0, x3,y3,0);

			neighbors->push_back(Neighbor(*itr, cost));
		}
	};

	// Recreate the start set.
	start_set.clear();
	for (int e = 0; e < mesh.number_of_edges(); ++e) {
		float y1 = mesh.get_point(mesh.get_edge(e).first).y;
		float y2 = mesh.get_point(mesh.get_edge(e).second).y;
		// Require that both y-coordinates are along the start line
		// because the first edge in the path is not counted.
		if (y1 < 0.5 || y2 < 0.5) {
			start_set.insert(e);
		}
	}

	std::cerr << "Edges and curvature\n";
	std::cerr << "----------------------\n";
	start_time = ::get_wtime();
	evaluations = 0;
	double dist = shortest_path(mesh.number_of_edges(), start_set, end_set, get_neighbors_curvature, &path);
	end_time = ::get_wtime();
	std::cerr << "Shortest path time: " << end_time - start_time << " seconds ";
	std::cerr << "(" << evaluations << " evaluations).\n";

	std::cerr << "Smallest curvature: " << dist << std::endl;
	std::cerr << "Path length: " << path.size() << std::endl << std::endl;

	mesh.draw_path(path, "0000ff");

	if (use_pairs) {
		// Compute the shortest path, but instead use edge pairs
		// as nodes in the graph.

		std::vector<int> adjacent_pairs;
		adjacent_pairs.reserve(1000);

		auto get_neighbors_pairs =
			[&mesh, &adjacent_pairs, &evaluations]
			(int ep, std::vector<Neighbor>* neighbors) -> void
		{
			evaluations++;

			mesh.get_adjacent_pairs(ep, &adjacent_pairs);
			for (auto itr = adjacent_pairs.begin(); itr != adjacent_pairs.end(); ++itr) {
				int p2, p3;
				tie(ignore, p2, p3) = mesh.get_edge_pair(*itr);

				float x2 = mesh.get_point(p2).x;
				float y2 = mesh.get_point(p2).y;
				float x3 = mesh.get_point(p3).x;
				float y3 = mesh.get_point(p3).y;

				float dx = x2 - x3;
				float dy = y2 - y3;

				float cost = std::sqrt(dx*dx + dy*dy);
				neighbors->push_back(Neighbor(*itr, cost));
			}
		};

		std::set<int> start_set_pairs, end_set_pairs;
		// Create the start and end set.
		for (int ep = 0; ep < mesh.number_of_edge_pairs(); ++ep) {
			int p1, p2, p3;
			tie(p1, p2, p3) = mesh.get_edge_pair(ep);

			float y1 = mesh.get_point(p1).y;
			float y2 = mesh.get_point(p2).y;
			float y3 = mesh.get_point(p3).y;

			// Require that all y-coordinates are on the start line
			// because the first pair in the path is not counted.
			if (y1 < 0.5 && y2 < 0.5 && y3 < 0.5) {
				start_set_pairs.insert(ep);
			}
			if (y1 > n - 1.5 || y2 > n - 1.5 || y3 > n - 1.5) {
				end_set_pairs.insert(ep);
			}
		}

		std::function<double(int)> heuristic_pairs =
			[n, &mesh, &heuristic_evaluations]
			(int i) -> double
		{
			heuristic_evaluations++;

			int p1, p2, p3;
			tie(p1, p2, p3) = mesh.get_edge_pair(i);

			float maxy = std::max(mesh.get_point(p1).y,
			             std::max(mesh.get_point(p2).y,
			                      mesh.get_point(p3).y));
			return std::max(0.0f, n - 1.5f - maxy);
		};
		//std::function<double(int)> heuristic_pairs(heuristic_lambda_pairs);

		// Perform the pair comparision
		std::vector<int> path_pairs;
		std::cerr << "Pairs\n";
		std::cerr << "-----------\n";
		compare_dijkstras_and_astar(mesh.number_of_edge_pairs(), get_neighbors_pairs,
									heuristic_pairs, start_set_pairs, end_set_pairs,
									&path_pairs, &evaluations, &heuristic_evaluations);

		// Convert pair path to edge path.
		path.clear();
		int e1, e2;
		for (auto itr = path_pairs.begin(); itr != path_pairs.end(); ++itr) {
			int p1, p2, p3;
			tie(p1, p2, p3) = mesh.get_edge_pair(*itr);

			e1 = mesh.find_edge(p1, p2);
			e2 = mesh.find_edge(p2, p3);
			path.push_back(e1);
		}
		path.push_back(e2);
		mesh.draw_path(path, "ffff00");
	}

	// Compute the shortest path, but instead use edge pairs
	// as nodes in the graph.

	auto get_neighbors_points =
		[&mesh, &evaluations]
		(int p, std::vector<Neighbor>* neighbors) -> void
	{
		evaluations++;

		auto& adjacent = mesh.get_point(p).adjacent_points;
		for (auto itr = adjacent.begin(); itr != adjacent.end(); ++itr) {
			float x1 = mesh.get_point(p).x;
			float y1 = mesh.get_point(p).y;
			float x2 = mesh.get_point(*itr).x;
			float y2 = mesh.get_point(*itr).y;

			float dx = x1 - x2;
			float dy = y1 - y2;

			float cost = std::sqrt(dx*dx + dy*dy);
			neighbors->push_back(Neighbor(*itr, cost));
		}
	};

	std::set<int> start_set_points, end_set_points;
	// Create the start and end set.
	for (int p = 0; p < mesh.number_of_points(); ++p) {
		float y = mesh.get_point(p).y;
		if (y < 0.05) {
			start_set_points.insert(p);
		}
		if (y > n - 1.05) {
			end_set_points.insert(p);
		}
	}

	std::function<double(int)> heuristic_points =
		[n, &mesh, &heuristic_evaluations]
		(int p) -> double
	{
		heuristic_evaluations++;

		float y = mesh.get_point(p).y;
		return y;
	};

	std::vector<int> path_points;
	std::cerr << "Points\n";
	std::cerr << "-----------\n";
	compare_dijkstras_and_astar(mesh.number_of_points(), get_neighbors_points,
								heuristic_points, start_set_points, end_set_points,
								&path_points, &evaluations, &heuristic_evaluations);
	return 0;
}

int main()
{
	try {
		return main_function();
	}
	catch (std::exception& e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		return 1;
	}
}
