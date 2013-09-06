// Petter Strandmark 2013.

#include <cmath>
#include <random>

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include <catch.hpp>
#include <curve_extraction/google_test_compatibility.h>

#include <curve_extraction/curvature.h>

using namespace curve_extraction;

TEST_CASE("compute_curvature/helix", "")
{
	double x1 = 1;
	double y1 = 0;
	double z1 = 0;

	double x2 = 0.885456025653210;
	double y2 = 0.464723172043769;
	double z2 = 0.025641025641026;

	double x3 = 0.568064746731156;
	double y3 = 0.822983865893656;
	double z3 = 0.051282051282051;

	// Real value   0.481288653714417
	// Approx value 0.495821512020759;

	const double power = 2.0;
	const int n_points = 1000000;
	double k2_int_pair = compute_curvature(x1,y1,z1,
	                                       x2,y2,z2,
	                                       x3,y3,z3,
	                                       power, false, n_points);
	CHECK(curvature_cache_misses == 1);
	EXPECT_NEAR(k2_int_pair, 0.495821512020759, 1e-6);

	// Iterate to test cache.
	for (int iter = 1; iter <= 10; ++iter) {
		k2_int_pair = compute_curvature(x1,y1,z1,
										x2,y2,z2,
										x3,y3,z3,
										power, true, n_points);
		CHECK(curvature_cache_misses == 2);
		CHECK(curvature_cache_hits == iter - 1);
	}
}

TEST_CASE("compute_torsion/helix", "")
{
	double x1 = 1;
	double y1 = 0;
	double z1 = 0;

	double x2 = 0.885456025653210;
	double y2 = 0.464723172043769;
	double z2 = 0.025641025641026;

	double x3 = 0.568064746731156;
	double y3 = 0.822983865893656;
	double z3 = 0.051282051282051;

	double x4 = 0.120536680255323;
	double y4 = 0.992708874098054;
	double z4 = 0.076923076923077;

	const int n_points = 1000000;

	// Iterate to test cache.
	for (int iter = 1; iter <= 10; ++iter) {
		double t1_int_pair = compute_torsion(x1,y1,z1,
											 x2,y2,z2,
											 x3,y3,z3,
											 x4,y4,z4,
											 1.0, n_points);
		EXPECT_NEAR(t1_int_pair, 0.026620592104080, 1e-6);
	}

	// Iterate to test cache.
	for (int iter = 1; iter <= 10; ++iter) {
		double t2_int_pair = compute_torsion(x1,y1,z1,
											 x2,y2,z2,
											 x3,y3,z3,
											 x4,y4,z4,
											 2.0, n_points);
		EXPECT_NEAR(t2_int_pair, 0.001522684996624, 1e-6);
	}
}

TEST_CASE("compute_torsion/crooked_line_in_plane", "")
{
	float test = compute_torsion<float>(3.0, 4.0, 4.0,
	                                    4.0, 3.0, 3.0,
	                                    3.0, 3.0, 1.0,
	                                    4.0, 2.0, 0.0);
	CHECK(test == test);
	EXPECT_NEAR(test, 0.0f, 1e-6f);
}
