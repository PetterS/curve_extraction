// Petter Strandmark 2013.
//
// Reconstructs a tree in 3D from a set of 2D images.
//
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
using std::ignore;
using std::tie;

#include <Eigen/Dense>

#define CHECK(arg) if (!(arg)) { throw std::runtime_error(#arg); }
#include <pgm_image.h>
using ceres::examples::PGMImage;

#include <curve_extraction/curvature.h>
#include <curve_extraction/data_term.h>
#include <curve_extraction/grid_mesh.h>
#include <curve_extraction/mesh.h>
#include <curve_extraction/shortest_path.h>

using namespace curve_extraction;
 
std::istream& operator>>(std::istream& in, Eigen::MatrixXd& M)
{
	for (int i = 0; i < M.rows(); ++i) {
		for (int j = 0; j < M.cols(); ++j) {
			std::cin >> M(i, j);
		}
	}
	return in;
}

std::istream& operator>>(std::istream& in, Eigen::Vector4d& X)
{
	for (int i = 0; i < 4; ++i) {
		std::cin >> X(i);
	}
	return in;
}

int main_function()
{
	using namespace std;

	cerr << "This program is best run with piped output and input, e.g.\n";
	cerr << "    tree_reconstruction < wintertree.txt > wintertree-out.txt\n";

	int number_of_images = -1;
	cin >> number_of_images;
	cerr << "Number of images: " << number_of_images << endl;

	// Original images.
	vector<PGMImage<double>> Is;
	// Data term images.
	vector<PGMImage<double>> Ds;
	// Projection matrices.
	vector<Eigen::MatrixXd> Ps;
	// 3D points.
	vector<Eigen::Vector4d> Xs;

	//
	// Read the input images and the data term images.
	//
	for (int i = 0; i < number_of_images; ++i) {
		string filename;
		cin >> filename;
		Is.push_back(PGMImage<double>(filename));
		if (Is.back().width() <= 0 || Is.back().height() <= 0) {
			throw runtime_error("Could not read image.");
		}
		cerr << filename << endl;

		cin >> filename;
		Ds.push_back(PGMImage<double>(filename));
		if (Ds.back().width() <= 0 || Ds.back().height() <= 0) {
			throw runtime_error("Could not read data image.");
		}

		Eigen::MatrixXd P(3, 4);
		cin >> P;
		cerr << filename << "\n" << P << "\n";
		Ps.push_back(P);
	}

	int number_of_points = -1;
	cin >> number_of_points;
	cerr << "Number of points: " << number_of_points << endl;

	//
	// Read the 3D points.
	//
	for (int i = 0; i < number_of_points; ++i) {
		Xs.push_back(Eigen::Vector4d());
		cin >> Xs[i];
		if (!cin) {
			throw runtime_error("Could not read point.");
		}
		Xs[i] /= Xs[i][3];
	}
	cerr << number_of_points << " points OK.\n";

	int root = 0;
	cin >> root;

	int number_of_endpoints = -1;
	cin >> number_of_endpoints;
	cerr << "Number of end points: " << number_of_endpoints << endl;
	vector<int> Es;
	for (int i = 0; i < number_of_endpoints; ++i) {
		int e;
		cin >> e;
		if (!cin) {
			throw runtime_error("Could not read point.");
		}
		Es.push_back(e);
	}

	if (false) {
		// Test to check coordinates.
		for (int i = 0; i < number_of_images; ++i) {
			// Project the root from 3D to the image
			// coordinates.
			Eigen::Vector3d x = Ps[i] * Xs[root];
			x[0] /= x[2];
			x[1] /= x[2];
			// Set the pixel to white.
			*Is[i].MutablePixel(x[0], x[1]) = 255;

			// Project the first end point.
			x = Ps[i] * Xs[Es[0]];
			x[0] /= x[2];
			x[1] /= x[2];
			// Set the pixel to white.
			*Is[i].MutablePixel(x[0], x[1]) = 255;

			// Write to temporary file.
			stringstream sout;
			sout << "t" << i + 1 << ".pgm";
			Is[i].WriteToFile(sout.str());
		}
	}

	Eigen::Vector3d min_point, max_point;
	for (int i = 0; i < 3; ++i) {
		min_point[i] = 1e100;
		max_point[i] = -1e100;
	}
	for (int i = 0; i < number_of_points; ++i) {
		min_point[0] = min(min_point[0], Xs[i][0]);
		min_point[1] = min(min_point[1], Xs[i][1]);
		min_point[2] = min(min_point[2], Xs[i][2]);

		max_point[0] = max(max_point[0], Xs[i][0]);
		max_point[1] = max(max_point[1], Xs[i][1]);
		max_point[2] = max(max_point[2], Xs[i][2]);
	}
	cerr << "Points range from (" << min_point.transpose() 
	     << ") to (" << max_point.transpose() << ")" << endl;

	//
	// The coordinate transformation from the mesh to 3D is
	// 
	//   x = offset[0] + resolution[0] * mesh_x.
	//   y = offset[1] + resolution[1] * mesh_y.
	//   z = offset[2] + resolution[2] * mesh_z.
	//
	int n = 50;
	for (int i = 0; i < 3; ++i) {
		min_point[i] -= 0.25;
		max_point[i] += 0.25;
	}
	Eigen::Vector3d offset = min_point;
	Eigen::Vector3d resolution = (max_point - min_point) / double(n - 1);
	// Square pixels in images.
	vector<double> voxel_dimensions(3, 1.0);

	// Create the mesh.
	GridMesh mesh(n, n, n, 4.0, false);
	cerr << "Mesh created with " << mesh.number_of_points() << " points and "
	     << mesh.number_of_edges() << " edges (" << mesh.get_connectivity()
	     << " connectivity)." << endl;

	// Transform the mesh coordinates.
	mesh.transform_points(offset[0], offset[1], offset[2],
	                      resolution[0], resolution[1], resolution[2]);
	cerr << "Mesh transformed." << endl;

	// This function converts an index to a point in the mesh to
	// a 2D point in an image. 
	// It requires the projection matrix and the image.
	auto mesh_index_to_image_point =
		[&mesh]
		(int p, const Eigen::MatrixXd& P, const PGMImage<double>& I)
		-> pair<double, double>
	{
		float mx = mesh.get_point(p).x;
		float my = mesh.get_point(p).y;
		float mz = mesh.get_point(p).z;
		Eigen::Vector3d x = P * Eigen::Vector4d(mx, my, mz, 1.0);
		x[0] /= x[2];
		x[1] /= x[2];

		// Check overflow.
		x[0] = max(0.0, x[0]);
		x[0] = min(double(I.width()-1), x[0]);
		x[1] = max(0.0, x[1]);
		x[1] = min(double(I.height()-1), x[1]);

		return make_pair(x[0], x[1]);
	};

	//
	// Neighborhood function for length regularization.
	//
	double length_regularization = 0;
	cin >> length_regularization;
	cerr << "Using length regularization " << length_regularization << endl;

	//
	// Data terms.
	//
	vector<std::unique_ptr<PieceWiseConstant>> data_terms;
	for (int i = 0; i < number_of_images; ++i) {
		data_terms.emplace_back(new PieceWiseConstant(&(Ds[i].data()[0]),
		                                           Ds[i].width(),
		                                           Ds[i].height(),
		                                           1,
		                                           voxel_dimensions));
	}

	auto get_neighbors_length =
		[&mesh,
		 number_of_images,
		 &Ps,
		 &Ds,
		 &data_terms,
		 &mesh_index_to_image_point,
		 &length_regularization
		]
		(int p, std::vector<Neighbor>* neighbors) -> void
	{
		const auto& adjacent = mesh.get_point(p).adjacent_points;
		for (auto itr = adjacent.begin(); itr != adjacent.end(); ++itr) {
			
			float dx = mesh.get_point(p).x - mesh.get_point(*itr).x;
			float dy = mesh.get_point(p).y - mesh.get_point(*itr).y;
			float dz = mesh.get_point(p).z - mesh.get_point(*itr).z;
			double length = sqrt(dx*dx + dy*dy + dz*dz);

			// Go through every image.
			double data_cost = 0.0;
			for (int i = 0; i < number_of_images; ++i) {
				auto coord1 = mesh_index_to_image_point(p, Ps[i], Ds[i]);
				auto coord2 = mesh_index_to_image_point(*itr, Ps[i], Ds[i]);

				double tmp_cost =
					data_terms[i]->evaluate_line_integral(
						coord1.first, coord1.second, 0.0,
				        coord2.first, coord2.second, 0.0);
				data_cost = max(tmp_cost, data_cost);
			}

			double cost = data_cost + length_regularization * length;

			neighbors->push_back(Neighbor(*itr, cost));
		}
	};

	// This function takes a 3D point and returns the index
	// of the closest point in the mesh.
	auto get_X_mesh_index =
		[&mesh, &offset, &resolution]
		(Eigen::Vector4d X)
		-> int
	{
		float x = (X[0] - offset[0]) / resolution[0];
		float y = (X[1] - offset[1]) / resolution[1];
		float z = (X[2] - offset[2]) / resolution[2];
		x = offset[0] + int(x + 0.5) * resolution[0];
		y = offset[1] + int(y + 0.5) * resolution[1];
		z = offset[2] + int(z + 0.5) * resolution[2];
		int int_X = mesh.find_point(x, y, z);
		if (int_X < 0) {
			throw runtime_error("Could not find closest point.");
		}
		return int_X;
	};

	// Plot mesh in each image.
	for (int i = 0; i < number_of_images; ++i) {
		// Make a copy of the image.
		auto I = Is[i];
		for (int p = 0; p < mesh.number_of_points(); ++p) {
			auto coord = mesh_index_to_image_point(p, Ps[i], Is[i]);
			try {
				*I.MutablePixel(coord.first, coord.second) = 255;
			}
			catch (...) {
			}
		}
		// Write to temporary file.
		stringstream sout;
		sout << "mesh" << i + 1 << ".pgm";
		
		if (!I.WriteToFile(sout.str())) {
			throw runtime_error("Could not write output mesh file.");
		}
		cerr << sout.str() << " written." << endl;
	}

	// Create output images.
	vector<PGMImage<double>> Iouts = Is;

	// This function takes a path and adds it to
	// the output images.
	auto add_path_to_output =
		[&Ps, &mesh_index_to_image_point]
		(const vector<int>& point_path,
		 vector<PGMImage<double>>& Is)
	{
		auto draw_pixel =
			[]
			(PGMImage<double>& I, int x, int y)
		{
			for (int dx = -2; dx <= 2; ++dx) {
				for (int dy = -2; dy <= 2; ++dy) {
					*I.MutablePixel(x + dx, y + dy) = 255;
				}
			}
		};

		// Add path to image.
		for (int i = 0; i < Is.size(); ++i) {
			for (int p_ind = 0; p_ind < point_path.size() - 1; ++p_ind) {
				auto coord1 = mesh_index_to_image_point(point_path[p_ind], Ps[i], Is[i]);
				auto coord2 = mesh_index_to_image_point(point_path[p_ind + 1], Ps[i], Is[i]);

				// Draw edge.
				for (double t = 0.0; t <= 1.0; t += 0.02) {
					double x = t * coord1.first + (1 - t) * coord2.first;
					double y = t * coord1.second + (1 - t) * coord2.second;
					*Is[i].MutablePixel(x, y) = 255;
				}

				// Draw start point.
				if (p_ind == 0) {
					draw_pixel(Is[i], coord1.first, coord1.second);
				}

				// Draw start point.
				if (p_ind == point_path.size() - 2) {
					draw_pixel(Is[i], coord2.first, coord2.second);
				}
			}
		}
	};

	ShortestPathOptions options;
	options.print_progress = true;

	// Set up the start set.
	set<int> start_set;
	int start_point = get_X_mesh_index(Xs[root]);
	start_set.insert(start_point);

	for (int end_point = 0; end_point < number_of_endpoints; ++end_point) {
		vector<int> path;

		set<int> end_set;
		int end_point_index = get_X_mesh_index(Xs[Es[end_point]]);
		end_set.insert(end_point_index);

		auto length_heuristic =
			[&mesh, &end_set, &length_regularization]
			(int p) -> double
		{
			double min_length = 1e100;

			for (auto itr = end_set.begin(); itr != end_set.end(); ++itr) {
				float dx = mesh.get_point(p).x - mesh.get_point(*itr).x;
				float dy = mesh.get_point(p).y - mesh.get_point(*itr).y;
				float dz = mesh.get_point(p).z - mesh.get_point(*itr).z;
				double length = sqrt(dx*dx + dy*dy + dz*dz);
				min_length = min(length, min_length);
			}
			return length_regularization * min_length;
		};


		cerr << "Computing shortest path for end point " << end_point + 1 << "...";
		shortest_path(mesh.number_of_points(), start_set, end_set, get_neighbors_length, &path, length_heuristic, options);
		cerr << endl;
		add_path_to_output(path, Iouts);

		// The path we found is part of the start set for the
		// next iteration.
		stringstream sout;
		sout << "w_" << end_point + 1 << ".path";
		ofstream fout(sout.str());
		for (int p : path) {
			start_set.insert(p);

			// Write points to standard output.
			cout << mesh.get_point(p).x << " " << mesh.get_point(p).y << " " << mesh.get_point(p).z << endl;
			// Write points to file.
			fout << mesh.get_point(p).x << " " << mesh.get_point(p).y << " " << mesh.get_point(p).z << endl;
		}
		cerr << sout.str() << " written." << endl;
	}

	// Write output images to disk.
	for (int i = 0; i < number_of_images; ++i) {
		stringstream sout;
		sout << "w" << length_regularization << "_" << i + 1 << ".pgm";
		if (!Iouts[i].WriteToFile(sout.str())) {
			throw runtime_error("Could not write output file.");
		}
		cerr << sout.str() << " written." << endl;
	}

	//
	// Neighborhood function for curvature regularization.
	//
	double curvature_regularization = 0;
	cin >> curvature_regularization;
	cerr << "Using curvature regularization " << curvature_regularization << endl;

	auto get_neighbors_curvature =
		[&mesh,
		 number_of_images,
		 &Ps,
		 &Ds,
		 &data_terms,
		 &mesh_index_to_image_point,
		 &length_regularization,
		 &curvature_regularization
		]
		(int e, std::vector<Neighbor>* neighbors) -> void
	{
		int p = mesh.get_edge(e).second;
		static std::vector<int> adjacent;
		mesh.get_adjacent_edges(e, &adjacent);
		for (auto itr = adjacent.begin(); itr != adjacent.end(); ++itr) {
			int p2 = mesh.get_edge(*itr).second;

			float dx = mesh.get_point(p).x - mesh.get_point(p2).x;
			float dy = mesh.get_point(p).y - mesh.get_point(p2).y;
			float dz = mesh.get_point(p).z - mesh.get_point(p2).z;
			double length = sqrt(dx*dx + dy*dy + dz*dz);

			// Go through every image.
			double data_cost = 0.0;
			for (int i = 0; i < number_of_images; ++i) {
				auto coord1 = mesh_index_to_image_point(p,  Ps[i], Ds[i]);
				auto coord2 = mesh_index_to_image_point(p2, Ps[i], Ds[i]);
				
				double tmp_cost =
					data_terms[i]->evaluate_line_integral(
						coord1.first, coord1.second, 0.0,
				        coord2.first, coord2.second, 0.0);
				data_cost = max(tmp_cost, data_cost);
			}

			int q1 = mesh.get_edge(e).first;
			int q2 = mesh.get_edge(e).second;
			int q3 = mesh.get_edge(*itr).second;
			float x1 = mesh.get_point(q1).x;
			float y1 = mesh.get_point(q1).y;
			float z1 = mesh.get_point(q1).z;
			float x2 = mesh.get_point(q2).x;
			float y2 = mesh.get_point(q2).y;
			float z2 = mesh.get_point(q2).z;
			float x3 = mesh.get_point(q3).x;
			float y3 = mesh.get_point(q3).y;
			float z3 = mesh.get_point(q3).z;
			float curvature = compute_curvature<float>(x1,y1,z1,
			                                           x2,y2,z2,
			                                           x3,y3,z3,
			                                           2.0);

			double cost = data_cost
			            + length_regularization * length
						+ curvature_regularization * curvature;

			neighbors->push_back(Neighbor(*itr, cost));
		}
	};

	if (curvature_regularization > 0) {
		// Create output images.
		vector<PGMImage<double>> ICouts = Is;

		// Set up the start set.
		set<int> start_set;
		int start_point = get_X_mesh_index(Xs[root]);
		for (int e = 0; e < mesh.number_of_edges(); ++e) {
			if (mesh.get_edge(e).first == start_point) {
				start_set.insert(e);
			}
		}

		for (int end_point = 0; end_point < number_of_endpoints; ++end_point) {
			vector<int> path;

			set<int> end_set;
			int end_point_index = get_X_mesh_index(Xs[Es[end_point]]);
			for (int e = 0; e < mesh.number_of_edges(); ++e) {
				if (mesh.get_edge(e).first == end_point_index) {
					end_set.insert(e);
				}
			}

			cerr << "Computing shortest path with curvature for end point " << end_point + 1 << "...";
			shortest_path(mesh.number_of_edges(), start_set, end_set, get_neighbors_curvature, &path, nullptr, options);
			cerr << endl;

			// The path we found is part of the start set for the
			// next iteration.
			std::vector<int> point_path;
			point_path.push_back(mesh.get_edge(path[0]).first);
			for (int e : path) {
				start_set.insert(e);
				point_path.push_back(mesh.get_edge(e).second);
			}

			stringstream sout;
			sout << "c_" << end_point + 1 << ".path";
			ofstream fout(sout.str());
			for (int p : point_path) {
				// Write points to standard output.
				cout << mesh.get_point(p).x << " " << mesh.get_point(p).y << " " << mesh.get_point(p).z << endl;
				// Write points to file.
				fout << mesh.get_point(p).x << " " << mesh.get_point(p).y << " " << mesh.get_point(p).z << endl;
			}
			cerr << sout.str() << " written." << endl;

			add_path_to_output(point_path, ICouts);
		}

		// Write output images to disk.
		for (int i = 0; i < number_of_images; ++i) {
			stringstream sout;
			sout << "c" << curvature_regularization << "_" << i + 1 << ".pgm";
			if (!ICouts[i].WriteToFile(sout.str())) {
				throw runtime_error("Could not write output file.");
			}
			cerr << sout.str() << " written." << endl;
		}
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
