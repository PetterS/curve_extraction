#pragma once 
// Note boundary goes in the middle of the pixels.
// In the PieceWiseConstant data cost the boundary is real voxel/pixel boundary.
typedef std::vector<std::pair<double,double> > Points_vector;
class Boundary_points
{
	public:
	  Boundary_points (const vector<double> voxel_dimensions, int M, int N)
	    :  voxel_dimensions(voxel_dimensions), M(M), N(N)
	 {}

	// This is similar to evaluate_line_integral 
	Points_vector operator () (double sx, double sy, double ex, double ey)
	{
		using std::abs;
		using std::sqrt;

		Points_vector points;

		double dx = (ex - sx)*voxel_dimensions[0];
		double dy = (ey - sy)*voxel_dimensions[1];
		double line_length = sqrt(dx*dx + dy*dy);

		typedef std::pair<double, int> crossing;

		// Integral invariant of direction
		auto intersection =
		[&]
		(int index_change, double line_length, double dx, double start,
		 std::vector<crossing>& crossings) -> void
		{
			double k = dx/line_length;

			if (abs(k) <= 1e-4f)
				return;

			auto shifted_start = start + 0.5;
			double current = 1.0 - fractional_part(shifted_start);

			if (k < 0) 
			{
				current = 1.0 - current;
				index_change = -index_change;
			}

			for (; current <= abs(dx); current = current + 1.0) 
				crossings.push_back( crossing(abs(current/dx), index_change) );
		};
	
		std::vector<crossing>* scratch_space;
		#ifdef USE_OPENMP
			// Maximum 100 threads.
			// Support for C++11 thread_local is not very good yet.
			static std::vector<crossing> local_space[100];
			scratch_space = &local_space[omp_get_thread_num()];
		#else
			static std::vector<crossing> local_space;
			scratch_space = &local_space;
		#endif

		scratch_space->clear();
		scratch_space->push_back( crossing(0.0, 0) );

		intersection(1,line_length, dx, sx, *scratch_space);
		intersection(M,line_length, dy, sy, *scratch_space);
		std::sort(scratch_space->begin(), scratch_space->end());

		points.push_back(std::pair<double,double>(sx,sy));

		// Remove small increments
		auto prev = scratch_space->begin();
		auto next = scratch_space->begin();
		next++;
		for (; next != scratch_space->end(); prev++, next++)
		{
			if ((next->first - prev->first) < 1e-5f)
				continue;

		  points.push_back( std::pair<double,double>(sx + (next->first)*dx, sy + (next->first)*dy) );  
		} 

		points.push_back(std::pair<double,double>(ex,ey));

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