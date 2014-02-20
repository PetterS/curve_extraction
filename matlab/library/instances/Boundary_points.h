#pragma once

// Note boundary goes in the middle of the pixels.
// In the PieceWiseConstant data cost the boundary is the voxel/pixel boundary.
typedef std::pair<double,double> boundary_point;

class Boundary_points
{
	public:
	  Boundary_points (const vector<double> voxel_dimensions, int M, int N)
	    :  voxel_dimensions(voxel_dimensions), M(M), N(N)
	 {}

	// This is similar to evaluate_line_integral
	std::vector<boundary_point> operator () (double sx, double sy, double ex, double ey)
	{
		using std::abs;
		using std::sqrt;

		double dx = (ex - sx)*voxel_dimensions[0];
		double dy = (ey - sy)*voxel_dimensions[1];
		double line_length = sqrt(dx*dx + dy*dy);

		// Invariant of direction
		auto intersection =
		[&]
		(int index_change, double line_length, double dx, double start,
		 std::vector<double>& crossings) -> void
		{
			double k = dx/line_length;

			if (abs(k) <= 1e-4f)
				return;

			double current;

			if (k > 0)
				current = 1.0 - fractional_part(start);
			else
				current = fractional_part(start);

			for (; current < abs(dx) - 1e-6; current = current + 1.0)
				crossings.push_back( abs(current/dx) );
		};

		std::vector<double>* scratch_space;
		#ifdef USE_OPENMP
			// Maximum 100 threads.
			// Support for C++11 thread_local is not very good yet.
			static std::vector<double> local_space[100];
			scratch_space = &local_space[omp_get_thread_num()];
		#else
			static std::vector<double> local_space;
			scratch_space = &local_space;
		#endif

		scratch_space->clear();
		scratch_space->push_back(0);
		intersection(1,line_length, dx, sx, *scratch_space);
		intersection(M,line_length, dy, sy, *scratch_space);
		std::sort(scratch_space->begin(), scratch_space->end());

		std::vector<boundary_point> points;
		points.push_back( boundary_point(sx,sy) );

		// Remove small increments
		auto prev = scratch_space->begin();
		auto next = scratch_space->begin();
		next++;
		for (; next != scratch_space->end(); prev++, next++)
		{
			if (abs(*next - *prev) < 1e-5)
				continue;

			points.push_back( boundary_point(sx + (*next)*dx,
																	sy + (*next)*dy)
											);
		}

		points.push_back( boundary_point(ex,ey) );

		return points;
	}

	double fractional_part(double x)
	{
		double integer_part;
		modf(x, &integer_part);
		return x - integer_part;
	}

	const vector<double> voxel_dimensions;
	int M; // rows
	int N; // cols
};
