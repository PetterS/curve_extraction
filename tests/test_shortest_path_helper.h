
#define TEST(function, type) TEST_CASE(#function "/" #type, "")

// I have problems with EXPECT_THROW and gcc.
// Seems to work with g++ 4.7.2.
#if 1 //|| defined(_MSC_VER) || defined(__clang__)
	TEST(shortest_path, invalid_start)
	{
		const int n = 100;
		auto get_neighbors =
			[n]
			(int i, std::vector<Neighbor>* neighbors) -> void
		{
			neighbors->push_back(Neighbor(n - 1, 1.0));
		};

		std::set<int> start_set;
		std::set<int> end_set;
		std::vector<int> path;
		end_set.insert(n - 1);
		EXPECT_THROW(shortest_path(n, start_set, end_set, get_neighbors, &path),
					 std::runtime_error);
		start_set.insert(-1);
		EXPECT_THROW(shortest_path(n, start_set, end_set, get_neighbors, &path),
					 std::runtime_error);

		start_set.clear();
		start_set.insert(0);
		EXPECT_NO_THROW(shortest_path(n, start_set, end_set, get_neighbors, &path));
	}

	TEST(shortest_path, invalid_end)
	{
		const int n = 100;
		auto get_neighbors =
			[n]
			(int i, std::vector<Neighbor>* neighbors) -> void
		{
			neighbors->push_back(Neighbor(n - 1, 1.0));
		};

		std::set<int> start_set;
		std::set<int> end_set;
		std::vector<int> path;
		start_set.insert(0);
		end_set.insert(-1);
		EXPECT_THROW(shortest_path(n, start_set, end_set, get_neighbors, &path),
					 std::runtime_error);
		end_set.clear();
		end_set.insert(n - 1);
		EXPECT_NO_THROW(shortest_path(n, start_set, end_set, get_neighbors, &path));
	}

	TEST(shortest_path, no_path)
	{
		auto get_neighbors =
			[]
			(int i, std::vector<Neighbor>* neighbors) -> void
		{
		};
		std::set<int> start_set;
		std::set<int> end_set;
		std::vector<int> path;
		start_set.insert(0);
		end_set.insert(9);
		EXPECT_THROW(shortest_path(10, start_set, end_set, get_neighbors, &path),
					 std::runtime_error);
	}

	TEST(shortest_path, negative_weights)
	{
		const int n = 10;
		double weight;
		auto get_neighbors =
			[&weight, n]
			(int i, std::vector<Neighbor>* neighbors) -> void
		{
			neighbors->push_back(Neighbor((i + 1) % n, weight));
		};
		std::set<int> start_set;
		std::set<int> end_set;
		std::vector<int> path;
		start_set.insert(0);
		end_set.insert(n - 1);

		weight = -1.0;
		EXPECT_THROW(shortest_path(n, start_set, end_set, get_neighbors, &path),
		             std::runtime_error);

		weight = 1.0;
		EXPECT_NO_THROW(shortest_path(n, start_set, end_set, get_neighbors, &path));
	}
#endif

TEST(shortest_path, line)
{
	const int n = 100;
	auto get_neighbors =
		[n]
		(int i, std::vector<Neighbor>* neighbors) -> void
	{
		neighbors->push_back(Neighbor((i + 1) % n, 1.0));
	};

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(0);
	end_set.insert(n - 1);
	double min_dist = shortest_path(n, start_set, end_set, get_neighbors, &path);
	EXPECT_NEAR(min_dist, double(n - 1), 1e-10);
	ASSERT_EQ(path.size(), n);
	for (int i = 0; i < n; ++i) {
		EXPECT_EQ(path[i], i);
	}
}

TEST(shortest_path, fully_connected)
{
	const int n = 100;
	auto get_neighbors =
		[n]
		(int i, std::vector<Neighbor>* neighbors) -> void
	{
		for (int j = 0; j < n; ++j) {
			if (i != j) {
				neighbors->push_back(Neighbor(j, 1.0));
			}
		}
	};

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(0);
	end_set.insert(n - 1);
	double min_dist = shortest_path(n, start_set, end_set, get_neighbors, &path);
	EXPECT_NEAR(min_dist, 1.0, 1e-10);
	ASSERT_EQ(path.size(), 2);
	EXPECT_EQ(path[0], 0);
	EXPECT_EQ(path[1], n - 1);
}

TEST(shortest_path, simple_grid4)
{
	//  0  1  2  3
	//  4  5  6  7
	//  8  9 10 11
	// 12 13 14 15
	// 16 17 18 19

	int n = 4;
	auto get_neighbors = [n](int i, std::vector<Neighbor>* neighbors) -> void {
		int x = i % n;
		int y = i / n;
		if (x > 0) {
			neighbors->push_back(Neighbor(i - 1, 1.0));
		}
		if (x < n - 1) {
			neighbors->push_back(Neighbor(i + 1, 1.0));
		}
		if (y > 0) {
			neighbors->push_back(Neighbor(i - n, 1.0));
		}
		if (y < n - 1) {
			neighbors->push_back(Neighbor(i + n, 1.0));
		}
	};

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(16);
	end_set.insert(11);
	double min_dist = shortest_path(20, start_set, end_set, get_neighbors, &path);
	EXPECT_NEAR(min_dist, 5.0, 1e-10);
	ASSERT_EQ(path.size(), 6);
	EXPECT_EQ(path[0], 16);
	EXPECT_EQ(path[5], 11);
}

TEST(shortest_path, simple_grid8)
{
	//  0  1  2  3
	//  4  5  6  7
	//  8  9 10 11
	// 12 13 14 15
	// 16 17 18 19

	int n = 4;
	auto get_neighbors = [n](int i, std::vector<Neighbor>* neighbors) -> void {
		int x = i % n;
		int y = i / n;
		if (x > 0) {
			neighbors->push_back(Neighbor(i - 1, 1.0));
		}
		if (x < n - 1) {
			neighbors->push_back(Neighbor(i + 1, 1.0));
		}
		if (y > 0) {
			neighbors->push_back(Neighbor(i - n, 1.0));
		}
		if (y < n - 1) {
			neighbors->push_back(Neighbor(i + n, 1.0));
		}

		if (x < n -1 && y < n - 1) {
			neighbors->push_back(Neighbor(i + 1 + n, std::sqrt(2.0)));
		}
		if (x > 0 && y < n - 1) {
			neighbors->push_back(Neighbor(i - 1 + n, std::sqrt(2.0)));
		}
		if (x < n -1 && y > 0) {
			neighbors->push_back(Neighbor(i + 1 - n, std::sqrt(2.0)));
		}
		if (x > 0 && y > 0) {
			neighbors->push_back(Neighbor(i - 1 - n, std::sqrt(2.0)));
		}
	};

	std::set<int> start_set;
	std::set<int> end_set;
	std::vector<int> path;
	start_set.insert(16);
	end_set.insert(10);
	end_set.insert(2);
	end_set.insert(3);
	end_set.insert(6);
	end_set.insert(7);
	double min_dist = shortest_path(20, start_set, end_set, get_neighbors, &path);
	EXPECT_NEAR(min_dist, 2.0 *  std::sqrt(2.0), 1e-6);
	ASSERT_EQ(path.size(), 3);
	EXPECT_EQ(path[0], 16);
	EXPECT_EQ(path[2], 10);
}
