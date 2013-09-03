// Petter Strandmark 2013.

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
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

#include <vessel/curvature.h>
#include <vessel/grid_mesh.h>
#include <vessel/shortest_path.h>

using namespace vessel;


int main_function()
{
	const int n = 10;
	bool use_pairs = true;

	// Coefficient in front of length.
	const double rho = 10000.0;

	// Coefficient in front of curvature.
	const double sigma = 10000.0;

	// Coefficient in front of torsion.
	const double nu = 10000.0;

	// Random seed.
	unsigned seed = 6;

	const double mesh_distance = 4.0;

	ShortestPathOptions options;
	options.print_progress = true;

	// This section defines the edges that are *not* part of the mesh.
	// All edges crossing a set of horizontal line segments are excluded.
	//
	// The line segments are  { y = y0
	//                        { x_min <= x <= x_max,
	// for y0 in {10.9, 10.1, 20.9, 20.1, 30.9, 30.1},
	//  x_min in {   0,    0, 2n/3,  2n/3,   0,    0} and
	//  x_max in { n/3,  n/3,    n,     n, n/3,  n/3}.
	std::vector<float> zvals;
	std::vector<float> xmin;
	std::vector<float> xmax;

	zvals.push_back(10.9f);
	xmin.push_back(2 * n / 3 + 0.1f);
	xmax.push_back(10000);
	zvals.push_back(10.1f);
	xmin.push_back(2 * n / 3 + 0.1f);
	xmax.push_back(10000);

	// This function returns whether a given edge should be
	// in the mesh.
	std::function<bool(float, float, float,
	                   float, float, float)>
		forbidden_edge = [n, &zvals, &xmin, &xmax]
		                  (float x1, float, float z1,
		                   float x2, float, float z2)
		               -> bool
	{
		for (int i = 0; i < zvals.size(); ++i) {
			float z0 = zvals[i];

			float dz = z2 - z1;
			if (std::abs(dz) > 0) {
				// z1 + t * dz = z0  <=>
				// t = (z0 - z1) / dz
				float t = (z0 - z1) / dz;
				if (0.0 < t && t < 1.0f) {
					float dx  = x2 - x1;
					float x0 = x1 + t * dx;
					if (x0 > xmax[i] && std::abs(dx) > 0.05) {
						return true;
					}
					if (x0 > xmax[i] + 0.2f) {
						return true;
					}
					if (x0 < xmin[i] && std::abs(dx) > 0.05) {
						return true;
					}
					if (x0 < xmin[i] - 0.2f) {
						return true;
					}
				}
			}
			else {
				if (std::abs(z1 - z0) < 0.05f &&
					x1 > xmax[i]) {
					return true;
				}
				if (std::abs(z2 - z0) < 0.05f &&
					x2 > xmax[i]) {
					return true;
				}
				if (std::abs(z1 - z0) < 0.05f &&
					x1 < xmin[i]) {
					return true;
				}
				if (std::abs(z2 - z0) < 0.05f &&
					x2 < xmin[i]) {
					return true;
				}
			}
		}
		return false;
	};

	// Add all allowed edges of length less than or
	// equal to 2.
	double start_time = ::get_wtime();

	// Create the mesh and add all n*n*n points.
	GridMesh mesh(n, n, n, mesh_distance, forbidden_edge, use_pairs);

	double end_time = ::get_wtime();
	std::cerr << "Create mesh time: " << end_time - start_time << " seconds.\n";

	std::cerr << mesh.number_of_points() << " points and "
	          << mesh.number_of_edges() << " edges.\n";

	std::cerr << double(mesh.number_of_edge_pairs()) / double(1e6) << " M edge pairs.\n\n";
	std::cerr << "Connectivity: " << mesh.get_connectivity() << std::endl;

	std::cerr << "Computing edge costs... ";
	std::vector<float> edge_cost(mesh.number_of_edges());
	for (int e = 0; e < mesh.number_of_edges(); ++e) {
		std::mt19937 engine((unsigned)e + seed);
		std::uniform_real_distribution<double> dist(0.0, 10.0);
		edge_cost[e] = dist(engine);
	}
	std::cerr << "done." << std::endl;

	// Main function for computing the shortest path. Given an edge e, it
	// fills the neighbors vector with a list of the neighbors and their
	// distances.
	int evaluations = 0;
	auto get_neighbors =
		[&mesh, &edge_cost, &evaluations, rho]
		(int e, std::vector<Neighbor>* neighbors) -> void
	{
		evaluations++;

		static std::vector<int> adjacent;
		mesh.get_adjacent_edges(e, &adjacent);
		for (auto itr = adjacent.begin(); itr != adjacent.end(); ++itr) {
			double cost = edge_cost[*itr] + rho * mesh.edge_length(*itr);
			neighbors->push_back(Neighbor(*itr, cost));
		}
	};

	// The two sets between which to compute the shortest path.
	std::set<int> start_set, end_set;
	// Create the start and end set.
	for (int e = 0; e < mesh.number_of_edges(); ++e) {
		float z1 = mesh.get_point(mesh.get_edge(e).first).z;
		float z2 = mesh.get_point(mesh.get_edge(e).second).z;
		// For length regularization, require that both
		// points lie in the plane.
		if (z1 < 0.5 && z2 < 0.5) {
			start_set.insert(e);
		}
		if (z1 > n - 1.5 || z2 > n - 1.5) {
			end_set.insert(e);
		}
	}

	// Compute the shortest path.
	evaluations = 0;
	std::vector<int> path;
	double dist;
	std::cerr << "Computing length... ";
	start_time = ::get_wtime();
	dist = shortest_path(mesh.number_of_edges(), start_set, end_set,
	                     get_neighbors, &path, 0, options);
	end_time = ::get_wtime();
	std::cerr << "done.\n";
	std::cerr << "Best cost : " << dist << " (" << path.size() << " elements in path)" << std::endl;
	std::cerr << "Length time: " << end_time - start_time << " seconds ";
	std::cerr << "(" << evaluations << " evaluations).\n";

	std::cerr << "Writing curve... ";
	std::ofstream fout("length.dat");
	if (path.size() > 0) {
		for (auto itr = path.begin(); itr != path.end(); ++itr) {
			const auto& e = mesh.get_edge(*itr);
			const auto& p2 = mesh.get_point(e.second);
			fout << p2.x << '\t' << p2.y << '\t' << p2.z << '\n';
		}
	}
	std::cerr << "done.\n\n";

	
	std::cerr << "Filling curvature cache... ";
	start_time = ::get_wtime();
	GridMesh cache_mesh(2 * mesh_distance, 2 * mesh_distance, 2 * mesh_distance, mesh_distance, true);

	float M = mesh_distance;
	int p1 = cache_mesh.find_point(M, M, M);
	const auto& adjacent_p = cache_mesh.get_point(p1).adjacent_points;

	for (auto itr2 = adjacent_p.begin(); itr2 != adjacent_p.end(); ++itr2) {
		int p2 = *itr2;
		for (auto itr3 = adjacent_p.begin(); itr3 != adjacent_p.end(); ++itr3) {
			int p3 = *itr3;

			float x = cache_mesh.get_point(p1).x;
			float y = cache_mesh.get_point(p1).y;
			float z = cache_mesh.get_point(p1).z;

			float x1 = 0;
			float y1 = 0;
			float z1 = 0;
			float x2 = cache_mesh.get_point(p2).x - x;
			float y2 = cache_mesh.get_point(p2).y - y;
			float z2 = cache_mesh.get_point(p2).z - z;
			float x3 = cache_mesh.get_point(p3).x - x + x2;
			float y3 = cache_mesh.get_point(p3).y - y + y2;
			float z3 = cache_mesh.get_point(p3).z - z + z2;
			compute_curvature<float>(x1,y1,z1, x2,y2,z2, x3,y3,z3, 2.0);
		}
	}
	end_time = ::get_wtime();
	std::cerr << "done in " << end_time - start_time << " seconds\n";
	std::cerr << "Curvature cache : "
	          << curvature_cache_hits/1000 << " khits / "
	          << curvature_cache_misses/1000 << " kmisses.\n\n"; 
	curvature_cache_hits = 0;
	curvature_cache_misses = 0;
	

	// Neighborhood function for curvature.
	auto get_neighbors_curvature =
		[&mesh, &edge_cost, &evaluations, sigma]
		(int e, std::vector<Neighbor>* neighbors) -> void
	{
		evaluations++;

		static std::vector<int> adjacent;
		mesh.get_adjacent_edges(e, &adjacent);

		neighbors->resize(adjacent.size());

		#ifdef USE_OPENMP
		#pragma omp parallel for
		#endif
		for (int i = 0; i < adjacent.size(); ++i) {
			int* itr = &adjacent[i];

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
			float z1 = mesh.get_point(p1).z;
			float x2 = mesh.get_point(p2).x;
			float y2 = mesh.get_point(p2).y;
			float z2 = mesh.get_point(p2).z;
			float x3 = mesh.get_point(p3).x;
			float y3 = mesh.get_point(p3).y;
			float z3 = mesh.get_point(p3).z;

			double cost = edge_cost[*itr] + sigma * compute_curvature<float>(x1,y1,z1, x2,y2,z2, x3,y3,z3, 2.0, false);
			(*neighbors)[i] = Neighbor(*itr, cost);
		}
	};

    start_set.clear();
    for (int e = 0; e < mesh.number_of_edges(); ++e) {
		float z1 = mesh.get_point(mesh.get_edge(e).first).z;
		float z2 = mesh.get_point(mesh.get_edge(e).second).z;
		// For curvature regularization, require that one
		// of the points lie in the plane.
		if (z1 < 0.5 || z2 < 0.5) {
			start_set.insert(e);
		}
		if (z1 > n - 1.5 || z2 > n - 1.5) {
			end_set.insert(e);
		}
	}

	// Compute the shortest path.
	evaluations = 0;
	path.clear();
	std::cerr << "Computing curvature... ";
	start_time = ::get_wtime();
	dist = shortest_path(mesh.number_of_edges(), start_set, end_set,
	                     get_neighbors_curvature, &path, 0, options);
	end_time = ::get_wtime();
	std::cerr << "done.\n";

	std::cerr << "Best cost : " << dist << " (" << path.size() << " elements in path)" << std::endl;
	std::cerr << "Curvature time: " << end_time - start_time << " seconds ";
	std::cerr << "(" << evaluations << " evaluations).\n";
	std::cerr << "Curvature cache : "
	          << curvature_cache_hits/1000 << " khits / "
	          << curvature_cache_misses/1000 << " kmisses.\n\n"; 

	std::cerr << "Writing curve... ";
	fout.close();
	fout.open("curvature.dat");
	if (path.size() > 0) {
		int p1;
		tie(p1, ignore) = mesh.get_edge(path[0]);
		const auto& point1 = mesh.get_point(p1);
		fout << point1.x << '\t' << point1.y << '\t' << point1.z << '\n';
		for (auto itr = path.begin(); itr != path.end(); ++itr) {
			int p2;
			tie(ignore, p2) = mesh.get_edge(*itr);
			const auto& point2 = mesh.get_point(p2);
			fout << point2.x << '\t' << point2.y << '\t' << point2.z << '\n';
		}
	}
	std::cerr << "done.\n\n";

	if (use_pairs) {
		//
		// Neighborhood function for torsion.
		//
		std::vector<int> adjacent_pairs;
		adjacent_pairs.reserve(10000);

		auto get_neighbors_torsion =
			[&mesh, &edge_cost, &evaluations, &adjacent_pairs, nu]
			(int ep, std::vector<Neighbor>* neighbors) -> void
		{
			evaluations++;

			mesh.get_adjacent_pairs(ep, &adjacent_pairs);

			for (auto itr = adjacent_pairs.begin(); itr != adjacent_pairs.end(); ++itr) {
				int q1, q2, q3, q4, q5, q6;
				tie(q1, q2, q3) = mesh.get_edge_pair(ep);
				tie(q4, q5, q6) = mesh.get_edge_pair(*itr);

				int p1, p2, p3, p4;
				if (q2 == q4 && q3 == q5) {
					p1 = q1;
					p2 = q2;
					p3 = q3;
					p4 = q6;
				}
				else {
					throw std::runtime_error("Mesh point error.");
				}

				int f1 = mesh.find_edge(q1, q2);
				int f2 = mesh.find_edge(q2, q3);
				int f3 = mesh.find_edge(q4, q5);
				int f4 = mesh.find_edge(q5, q6);
				int e3;
				if (f2 == f3) {
					e3 = f4;
				}
				else {
					throw std::runtime_error("Mesh pair error");
				}

				float x1 = mesh.get_point(p1).x;
				float y1 = mesh.get_point(p1).y;
				float z1 = mesh.get_point(p1).z;
				float x2 = mesh.get_point(p2).x;
				float y2 = mesh.get_point(p2).y;
				float z2 = mesh.get_point(p2).z;
				float x3 = mesh.get_point(p3).x;
				float y3 = mesh.get_point(p3).y;
				float z3 = mesh.get_point(p3).z;
				float x4 = mesh.get_point(p4).x;
				float y4 = mesh.get_point(p4).y;
				float z4 = mesh.get_point(p4).z;

				double cost = edge_cost[e3] + nu *
					compute_torsion<float>(
						x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, 2.0, 50);
				neighbors->push_back(Neighbor(*itr, cost));
			}
		};

		std::cerr << "Creating start/end sets for torsion... ";
		// The two sets between which to compute the shortest path.
		std::set<int> start_set_pairs, end_set_pairs;
		// Create the start and end set.
		for (int ep = 0; ep < mesh.number_of_edge_pairs(); ++ep) {
			int q1, q2, q3;
			tie(q1, q2, q3) = mesh.get_edge_pair(ep);

			int e1 = mesh.find_edge(q1 ,q2);
			int e2 = mesh.find_edge(q2, q3);
			if (start_set.count(e1) > 0 ||
				start_set.count(e2) > 0) {
					start_set_pairs.insert(ep);
			}
			if (end_set.count(e1) > 0 ||
				end_set.count(e2) > 0) {
					end_set_pairs.insert(ep);
			}
		}
		std::cerr << "done.\n";

		// Compute the shortest path with torsion.
		evaluations = 0;
		path.clear();
		std::cerr << "Computing torsion... ";
		start_time = ::get_wtime();
		dist = shortest_path(mesh.number_of_edge_pairs(), start_set_pairs, end_set_pairs,
							 get_neighbors_torsion, &path, 0, options);
		end_time = ::get_wtime();
		std::cerr << "done.\n";

		std::cerr << "Best cost     : " << dist << " (" << path.size() << " elements in path)" << std::endl;
		std::cerr << "Torsion time  : " << end_time - start_time << " seconds ";
		std::cerr << "(" << evaluations << " evaluations).\n";
		std::cerr << "Torsion cache : " 
				  << torsion_cache_hits / 1000 << " / "
				  << torsion_cache_misses / 1000
				  << " khits / kmisses.\n"; 

		std::cerr << "Writing curve... ";
		fout.close();
		fout.open("torsion.dat");
		if (path.size() > 0) {
			int p1, p2;
			tie(p1, p2, ignore) = mesh.get_edge_pair(path[0]);
			const auto& point1 = mesh.get_point(p1);
			const auto& point2 = mesh.get_point(p2);
			fout << point1.x << '\t' << point1.y << '\t' << point1.z << '\n';
			fout << point2.x << '\t' << point2.y << '\t' << point2.z << '\n';
			for (auto itr = path.begin(); itr != path.end(); ++itr) {
				int p3;
				tie(ignore, ignore, p3) = mesh.get_edge_pair(*itr);
				const auto& point3 = mesh.get_point(p3);
				fout << point3.x << '\t' << point3.y << '\t' << point3.z << '\n';
			}
		}
		std::cerr << "done.\n";
	}

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
