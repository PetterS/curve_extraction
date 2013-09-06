// Petter Strandmark 2013.

#include <random>
#include <stdexcept>

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include <catch.hpp>
#include <curve_extraction/google_test_compatibility.h>


#include <curve_extraction/shortest_path.h>

using namespace curve_extraction;

// First, include the tests for the "shortest_path" function.
#include "test_shortest_path_helper.h"

#if 0
// Then, include the tests for the "shortest_path_memory_efficient"
// function.
#define shortest_path shortest_path_memory_efficient
#include "test_shortest_path_helper.h"
#undef shortest_path
#endif

TEST_CASE("shortest_path/maximum_queue_size", "")
{
	//  0  1  2  3
	//  4  5  6  7
	//  8  9 10 11
	// 12 13 14 15
	// 16 17 18 19

	int m = 4;
	int n = 5;
	auto get_neighbors = [m, n](int i, std::vector<Neighbor>* neighbors) -> void {
		int x = i % m;
		int y = i / m;
		if (x > 0) {
			neighbors->push_back(Neighbor(i - 1, 1.0));
		}
		if (x < m - 1) {
			neighbors->push_back(Neighbor(i + 1, 1.0));
		}
		if (y > 0) {
			neighbors->push_back(Neighbor(i - m, 1.0));
		}
		if (y < n - 1) {
			neighbors->push_back(Neighbor(i + m, 1.0));
		}
	};

	ShortestPathOptions options;
	options.maximum_queue_size = 2;

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(16);
	end_set.insert(3);
	EXPECT_THROW(
		double min_dist = shortest_path(m * n, start_set, end_set, get_neighbors, &path, 0, options),
		std::runtime_error);
}


TEST_CASE("shortest_path/all_distances", "")
{
	//  0  1  2  3
	//  4  5  6  7
	//  8  9 10 11
	// 12 13 14 15
	// 16 17 18 19

	int m = 4;
	int n = 5;
	auto get_neighbors = [m, n](int i, std::vector<Neighbor>* neighbors) -> void {
		int x = i % m;
		int y = i / m;
		if (x > 0) {
			neighbors->push_back(Neighbor(i - 1, 1.0));
		}
		if (x < m - 1) {
			neighbors->push_back(Neighbor(i + 1, 1.0));
		}
		if (y > 0) {
			neighbors->push_back(Neighbor(i - m, 1.0));
		}
		if (y < n - 1) {
			neighbors->push_back(Neighbor(i + m, 1.0));
		}
	};

	ShortestPathOptions options;
	options.compute_all_distances = true;

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(13);
	shortest_path(m * n, start_set, end_set, get_neighbors, &path, 0, options);

	ASSERT_EQ(options.distance.size(), m * n);
	EXPECT_FLOAT_EQ(options.distance[0], 4.0f);
	EXPECT_FLOAT_EQ(options.distance[1], 3.0f);
	EXPECT_FLOAT_EQ(options.distance[2], 4.0f);
	EXPECT_FLOAT_EQ(options.distance[3], 5.0f);
	EXPECT_FLOAT_EQ(options.distance[4], 3.0f);
	EXPECT_FLOAT_EQ(options.distance[5], 2.0f);
	EXPECT_FLOAT_EQ(options.distance[6], 3.0f);
	EXPECT_FLOAT_EQ(options.distance[7], 4.0f);
	EXPECT_FLOAT_EQ(options.distance[8], 2.0f);
	EXPECT_FLOAT_EQ(options.distance[9], 1.0f);

	EXPECT_FLOAT_EQ(options.distance[10], 2.0f);
	EXPECT_FLOAT_EQ(options.distance[11], 3.0f);
	EXPECT_FLOAT_EQ(options.distance[12], 1.0f);
	EXPECT_FLOAT_EQ(options.distance[13], 0.0f);
	EXPECT_FLOAT_EQ(options.distance[14], 1.0f);
	EXPECT_FLOAT_EQ(options.distance[15], 2.0f);
	EXPECT_FLOAT_EQ(options.distance[16], 2.0f);
	EXPECT_FLOAT_EQ(options.distance[17], 1.0f);
	EXPECT_FLOAT_EQ(options.distance[18], 2.0f);
	EXPECT_FLOAT_EQ(options.distance[19], 3.0f);
}

TEST_CASE("A_star/perfect_heuristic", "")
{
	//  0  1  2  3
	//  4  5  6  7
	//  8  9 10 11
	// 12 13 14 15
	// 16 17 18 19

	int m = 4;
	int n = 5;
	int evaluations = 0;
	auto get_neighbors =
		[m, n, &evaluations]
		(int i, std::vector<Neighbor>* neighbors) -> void
	{
		evaluations++;

		int x = i % m;
		int y = i / m;
		if (x > 0) {
			neighbors->push_back(Neighbor(i - 1, 1.0));
		}
		if (x < m - 1) {
			neighbors->push_back(Neighbor(i + 1, 1.0));
		}
		if (y > 0) {
			neighbors->push_back(Neighbor(i - m, 1.0));
		}
		if (y < n - 1) {
			neighbors->push_back(Neighbor(i + m, 1.0));
		}
	};

	// Function which uses a perfect heuristic. This will
	// make A* go straight to the goal.
	auto heuristic_lambda =
		[m, n]
		(int i) -> double
	{
		int x = i % m;
		int y = i / m;
		return std::abs(y - 4) + std::abs(x - 3);
	};
	std::function<double(int)> heuristic(heuristic_lambda);

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(16);
	end_set.insert(19);

	evaluations = 0;
	shortest_path(m * n, start_set, end_set, get_neighbors, &path);
	// Regular Dijkstra's reduces to BFS in this case.
	EXPECT_GE(evaluations, 9);

	evaluations = 0;
	double min_dist = shortest_path(20, start_set, end_set,
	                                get_neighbors, &path,
	                                &heuristic);
	EXPECT_EQ(evaluations, 3);
	EXPECT_NEAR(min_dist, 3.0, 1e-10);
	ASSERT_EQ(path.size(), 4);
	EXPECT_EQ(path[0], 16);
	EXPECT_EQ(path[3], 19);
}

TEST_CASE("A_star/random_grid", "")
{
	const int n = 100;
	const double min_weight = 1.0;
	const double max_weight = 2.0;

	int evaluations = 0;
	auto get_neighbors =
		[n, &evaluations, min_weight, max_weight]
		(int i, std::vector<Neighbor>* neighbors) -> void
	{
		evaluations++;

		std::mt19937 engine((unsigned)i);
		std::uniform_real_distribution<double> rand_dist(min_weight, max_weight);
		auto rand = std::bind(rand_dist, engine);

		int x = i % n;
		int y = i / n;
		if (x > 0) {
			neighbors->push_back(Neighbor(i - 1, rand()));
		}
		if (x < n - 1) {
			neighbors->push_back(Neighbor(i + 1, rand()));
		}
		if (y > 0) {
			neighbors->push_back(Neighbor(i - n, rand()));
		}
		if (y < n - 1) {
			neighbors->push_back(Neighbor(i + n, rand()));
		}
		if (x < n - 1 && y < n - 1) {
			neighbors->push_back(Neighbor(i + 1 + n, rand()));
		}
	};

	// Heuristic assuming all edges have min_weight. This will
	// obviously not overestimate the shortest path length.
	auto heuristic_lambda =
		[n, min_weight]
		(int i) -> double
	{
		int x = i % n;
		int y = i / n;
		return min_weight * std::max(std::abs(y - (n - 1)), std::abs(x - (n - 1)));
	};
	std::function<double(int)> heuristic(heuristic_lambda);

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(0);
	end_set.insert(n*n - 1);

	evaluations = 0;
	double min_dist = shortest_path(n*n, start_set, end_set, get_neighbors, &path);
	int evaluations_dijkstra = evaluations;

	evaluations = 0;
	double min_dist_A_star = shortest_path(n*n, start_set, end_set,
	                                       get_neighbors, &path,
	                                       heuristic);
	EXPECT_LT(evaluations, evaluations_dijkstra);
	EXPECT_LT( std::abs(min_dist - min_dist_A_star) / min_dist, 1e-6);

	// Check that the path sums correctly.
	double d = 0;
	for (int i = 0; i < path.size() - 1; ++i) {
		// Get the neighbors.
		std::vector<Neighbor> neighbors;
		get_neighbors(path[i], &neighbors);
		// Find the next node along the path.
		auto itr = neighbors.begin();
		for (; itr != neighbors.end(); ++itr) {
			if (itr->destination == path[i + 1]) {
				break;
			}
		}
		// Make sure it was found.
		ASSERT_LT(itr, neighbors.end());
		// Add the distance along the edge.
		d += itr->distance;
	}
	EXPECT_LT( std::abs(min_dist - d) / min_dist, 1e-6);
}
