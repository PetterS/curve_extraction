// Petter Strandmark 2013.
//
// Outputs an SVG file "curve_completion.svg"
//
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <tuple>
using std::ignore;
using std::tie;


#include <vessel/curvature.h>
#include <vessel/grid_mesh.h>
#include <vessel/mesh.h>
#include <vessel/shortest_path.h>

using namespace vessel;


int main_function()
{
	const int n = 30;

	// Create the mesh
	GridMesh mesh(n, n, 6.0);

	double length_coefficient = 1.0;
	double curvature_coefficient = 1.0;
	double curvature_power = 1.0;

	// Neighborhood function for curvature computation.
	auto get_neighbors_curvature =
		[&mesh, &curvature_power, &curvature_coefficient, &length_coefficient]
		(int e, std::vector<Neighbor>* neighbors) -> void
	{
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

			float cost = curvature_coefficient 
			             * compute_curvature<float>(x1,y1,0,
			                                        x2,y2,0,
			                                        x3,y3,0,
			                                        curvature_power);
			cost += length_coefficient * (std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
				                          + std::sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2)));

			neighbors->push_back(Neighbor(*itr, cost));
		}
	};

	std::set<int> start_set, end_set;
	start_set.insert(mesh.find_edge(mesh.find_point(3*n/4, 0, 0),
	                                mesh.find_point(3*n/4, 1, 0)));
	end_set.insert(mesh.find_edge(mesh.find_point(1, 3*n/4, 0),
	                              mesh.find_point(0, 3*n/4, 0)));

	std::vector<int> path;
	ShortestPathOptions options;
	options.print_progress = true;
	mesh.start_SVG("curve_completion.svg",
	               [](int){return "stroke-width:0.005;stroke:#4f4f4f;";});



	length_coefficient    = 0.0;
	curvature_coefficient = 1.0;
	curvature_power = 1.0;
	std::cerr << "Computing power = " << curvature_power;
	shortest_path(mesh.number_of_edges(), start_set, end_set, get_neighbors_curvature, &path, 0, options);
	std::cerr << std::endl;
	mesh.draw_path(path, "ff0000");

	length_coefficient    = 0.0;
	curvature_coefficient = 1.0;
	curvature_power = 2.0;
	std::cerr << "Computing power = " << curvature_power;
	shortest_path(mesh.number_of_edges(), start_set, end_set, get_neighbors_curvature, &path, 0, options);
	std::cerr << std::endl;
	mesh.draw_path(path, "00ff00");

	length_coefficient    = 1.0;
	curvature_coefficient = 0.0;
	curvature_power       = 0.0;
	std::cerr << "Computing power = " << curvature_power;
	shortest_path(mesh.number_of_edges(), start_set, end_set, get_neighbors_curvature, &path, 0, options);
	std::cerr << std::endl;
	mesh.draw_path(path, "0000ff");

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
