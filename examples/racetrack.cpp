// Petter Strandmark 2014.
//
// Racetrack game solver.
//
// See https://en.wikipedia.org/wiki/Racetrack_(game)

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <tuple>
using namespace std;

#include <curve_extraction/shortest_path.h>
using namespace curve_extraction;

#define CHECK(arg) if (!(arg)) { throw std::runtime_error(#arg); }
#include <pgm_image.h>
using ceres::examples::PGMImage;

void print_field(
	const vector<string>& field,
	const vector<pair<int, int>>& path)
{
	const int griddim = 30;
	PGMImage<int> image(30 * field.at(0).length(), 30 * field.size());
	image.Set(255);

	for (size_t i = 0; i < field.size(); ++i) {
	for (size_t j = 0; j < field.at(0).length(); ++j) {
		int x = griddim / 2 + griddim * j;
		int y = griddim / 2 + griddim * i;

		if (field.at(i).at(j) != 'X') {
			for (int ind = -griddim / 2; ind < griddim / 2; ++ind) {
				*image.MutablePixel(x + ind, y) = 200;
				*image.MutablePixel(x, y + ind) = 200;
			}
		}
	}}

	for (size_t ind = 0; ind < path.size(); ++ind) {
		auto& coord = path[ind];
		int x = griddim / 2 + griddim * coord.second;
		int y = griddim / 2 + griddim * coord.first;
		for (int dx = -2; dx <= +2; ++dx) {
		for (int dy = -2; dy <= +2; ++dy) {
			*image.MutablePixel(x+dx, y+dy) = 0;
		}}

		if (ind > 0) {
			auto& prev_coord = path[ind - 1];
			int x0 = griddim / 2 + griddim * prev_coord.second;
			int y0 = griddim / 2 + griddim * prev_coord.first;
			for (double t = 0.0; t <= 1.0; t+=0.01) {
				int xp = x + (x0 - x) * t + 0.5;
				int yp = y + (y0 - y) * t + 0.5;
				*image.MutablePixel(xp, yp) = 0;
			}
		}
	}

	clog << "Writing file to field.pgm." << endl;
	image.WriteToFile("field.pgm");
}

int main_function()
{
	//
	// This vector defines the playing board. 
	// 
	//  'E' is the end set.
	//  'S' is the start set.
	//  'X' can not be passed.
	//
	vector<string> field1;
	field1.emplace_back("XXXXXXXXXXXXXXXXXXXXX");
	field1.emplace_back("X             XXXXXXX");
	field1.emplace_back("X XXXXXX          XXX");
	field1.emplace_back("X     XXXXX        XX");
	field1.emplace_back("X       XXXXXXXX   XX");
	field1.emplace_back("X         XXXXXX    X");
	field1.emplace_back("X         X         X");
	field1.emplace_back("X        EXS        X");
	field1.emplace_back("XXXXXXXXXXXXXXXXXXXXX");

	// Another level.
	vector<string> field2;
	field2.emplace_back("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
	field2.emplace_back("XX                            XX       XX");
	field2.emplace_back("XX  XXXXXXXXXXXXXXXXXXXXXXXX   X       XX");
	field2.emplace_back("XXX     XXXXXXXXXXXXXXXX      XX   XX  XX");
	field2.emplace_back("XXXX         XXXXXXXXXX      XXX  XXX  XX");
	field2.emplace_back("XXXXXX           XXXXXX     XXXX  XXX  XX");
	field2.emplace_back("XXXXXXXXXX        XXXXX   XXXXXX  XXX  XX");
	field2.emplace_back("XXXXXXXXXXXX       XXX     XXXXX  XXX  XX");
	field2.emplace_back("XXXXXXXXXXXXXXXX   XXXXX    XXXX  XXX  XX");
	field2.emplace_back("XXXXXXXXXX          XXXX     XXX  XXX  XX");
	field2.emplace_back("XX                 XXXXXX     XX  XXX  XX");
	field2.emplace_back("XX               XXXXXXXX      X  XXX  XX");
	field2.emplace_back("XX  XXXXXXXXXXXXXXXXXXXXXX     X  XXX  XX");
	field2.emplace_back("XX  XXXXXXXXXXXXXXXXXXXXXXXX     XXXX  XX");
	field2.emplace_back("XX   XX    XXXXXX  EXS  XXXXXXXXXXXX   XX");
	field2.emplace_back("XXX   X            EXS    XXXXXXXXXX   XX");
	field2.emplace_back("XXXX        XXX    EXS                 XX");
	field2.emplace_back("XXXXXX    XXXXXXX  EXS                 XX");
	field2.emplace_back("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");

	// Choose level.
	const auto& field = field2;

	auto M = field.size();
	auto N = field.at(0).size();

	clog << "Track has dimensions " << M << "x" << N << endl;

	auto coord_to_line_graph_node = [N, M](int i1, int j1, int i2, int j2) -> int
	{
		if (i1 < 0 || i1 >= M) return -1;
		if (i2 < 0 || i2 >= M) return -1;
		if (j1 < 0 || j1 >= N) return -1;
		if (j2 < 0 || j2 >= N) return -1;

		int ind1 = i1 + M * j1;
		int ind2 = i2 + M * j2;
		return ind1 + N * M * ind2;
	};

	auto line_graph_node_to_coord = [N, M](int ind)
		-> tuple<int, int, int, int>
	{
		int ind1 = ind % (N * M);
		int ind2 = ind / (N * M);
		int i1 = ind1 % M;
		int j1 = ind1 / M;
		int i2 = ind2 % M;
		int j2 = ind2 / M;
		return make_tuple(i1, j1, i2, j2);
	};

	auto valid_edge = [N, M, &field](int i1, int j1, int i2, int j2)
		-> bool
	{
		if (field.at(i1).at(j1) == 'X') return false;
		if (field.at(i2).at(j2) == 'X') return false;
		for (double t = 0.0; t <= 1.0; t+=0.01) {
			int i = i1 + (i2 - i1) * t + 0.5;
			int j = j1 + (j2 - j1) * t + 0.5;
			if (field.at(i).at(j) == 'X') return false;
		}
		return true;
	};

	set<int> start_set;
	for (size_t i = 0; i < M; ++i) {
	for (size_t j = 0; j < N; ++j) {
		if (field.at(i).at(j) == 'S') {
			for (int di = -1; di <= +1; ++di) {
			for (int dj = -1; dj <= +1; ++dj) {
			if (di != 0 || dj != 0) {
				if (field.at(i + di).at(j + dj) != 'X') {
					start_set.emplace(coord_to_line_graph_node(i, j, i+di, j+dj));
				}
			}}}
		}
	}}

	set<int> end_set;
	for (size_t i = 0; i < M; ++i) {
	for (size_t j = 0; j < N; ++j) {
		if (field.at(i).at(j) == 'E') {
			for (int di = -1; di <= +1; ++di) {
			for (int dj = -1; dj <= +1; ++dj) {
			if (di != 0 || dj != 0) {
				if (field.at(i + di).at(j + dj) != 'X') {
					end_set.emplace(coord_to_line_graph_node(i+di, j+dj, i, j));
				}
			}}}
		}
	}}

	auto get_neighbors =
		[M, 
		 N,
		 &field,
		 &line_graph_node_to_coord,
		 &coord_to_line_graph_node,
		 &valid_edge]
		(int e, vector<Neighbor>* neighbors) -> void
	{
		auto coords = line_graph_node_to_coord(e);
		auto from_i = get<0>(coords);
		auto from_j = get<1>(coords);
		auto to_i = get<2>(coords);
		auto to_j = get<3>(coords);

		auto i = to_i + (to_i - from_i);
		auto j = to_j + (to_j - from_j);

		for (int di = -1; di <= +1; ++di) {
		for (int dj = -1; dj <= +1; ++dj) {
			int i2 = i+di;
			int j2 = j+dj;
			if (i2 > 0 && i2 < M && j2 > 0 && j2 < N &&
				field.at(i2).at(j2) != 'X' &&
				valid_edge(to_i, to_j, i2, j2)) {
				neighbors->emplace_back(coord_to_line_graph_node(to_i, to_j, i2, j2), 1.0);
			}
		}}
	};

	vector<int> path;
	shortest_path(M*N*M*N, start_set, end_set, get_neighbors, &path);
	clog << path.size() << " elements in path." << endl;

	vector<pair<int, int>> point_path;
	auto first_coords = line_graph_node_to_coord(path.at(0));
	point_path.emplace_back(get<0>(first_coords), get<1>(first_coords));
	for (int ind: path) {
		auto coords = line_graph_node_to_coord(ind);
		point_path.emplace_back(get<2>(coords), get<3>(coords));
	}

	print_field(field, point_path);

	return 0;
}

int main()
{
	try {
		return main_function();
	}
	catch (exception& e) {
		cerr << "ERROR: " << e.what() << endl;
		return 1;
	}
}
