#pragma once
#include <spii-thirdparty/fadiff.h>
#include <spii/auto_diff_term.h>
namespace fadbad
{

template<unsigned n>
F<double, n> abs(F<double, n> x)
{
	if (spii::to_double(x) < 0) {
		return -x;
	}
	else {
		return x;
	}
}

template<unsigned n>
F<F<double, n>, n> abs(F<F<double, n>, n> x)
{
	if (spii::to_double(x) < 0) {
		return -x;
	}
	else {
		return x;
	}
}
}

template<typename R> 
R fractional_part(R x) 
{
	double integer_part;
	modf(spii::to_double(x), &integer_part);
	return x - integer_part;
}

// Container
template<typename R>
class Boundary_points
{
public:
	Boundary_points(const matrix<double>& data, const vector<double>& vd, R _x, R _y) 
	: data(data), vd(vd)
	{
		x.push_back(_x);
		y.push_back(_y);
	};

	int num_triplets()
	{
		return num_pairs()-1;
	}

	int num_pairs()
	{
		return linear_indices.size();
	}

	typedef std::tuple<R,R,R,R,R,R, double, double, double> triplet_tuple;
	triplet_tuple get_triplet(int triplet)
	{	
		int y_int = linear_indices[triplet + 1]/data.M;
		int x_int = linear_indices[triplet + 1] - y_int*data.M;

		double d11,d10,d01;
		tie(d11,d10,d01) = get_coefficents(x_int,y_int);

		R x1 = x[triplet + 0];
		R x2 = x[triplet + 1];
		R x3 = x[triplet + 2];
		R y1 = y[triplet + 0];
		R y2 = y[triplet + 1];
		R y3 = y[triplet + 2]; 

		return triplet_tuple(x1,x2,x3,y1,y2,y3,d11,d10,d01);
	}

	typedef std::tuple<R,R,R,R, double, double, double> pair_tuple;
	pair_tuple get_pair(int pair)
	{
		int y_int = linear_indices[pair]/data.M;
		int x_int = linear_indices[pair] - y_int*data.M;

		double d11,d10,d01;
		tie(d11,d10,d01) = get_coefficents(x_int,y_int);

		R dx = (x[pair+1] - x[pair])*vd[0];
    R dy = (y[pair+1] - y[pair])*vd[1];

	  R x0 = (x[pair] - x_int)	*vd[0];
  	R y0 = (y[pair] - y_int)	*vd[1];

    return pair_tuple(x0,y0,dx,dy, d11, d10, d01);
	}

	typedef std::tuple<double, double, double> coefficents;
	coefficents get_coefficents(int x_int, int y_int)
	{
		double i00 = data_value(x_int + 0, y_int + 0)*vd[2];
		double i01 = data_value(x_int + 0, y_int + 1)*vd[2];
		double i10 = data_value(x_int + 1, y_int + 0)*vd[2];
		double i11 = data_value(x_int + 1, y_int + 1)*vd[2];

		double d11 = (i00 + i11 - i10- i01)/(vd[0]*vd[1]);
		double d10 = (i10 - i00)/(vd[0]*vd[1]);
		double d01 = (i01 - i00)/(vd[0]*vd[1]);

		return coefficents(d11,d10,d01);
	}

	// All subsequent points
	void add(R _x, R _y, int linear_index) 
	{
		x.push_back(_x);
		y.push_back(_y);
		linear_indices.push_back(linear_index);
  }

  // Sample data or closest point inside data.
  double data_value(int x, int y) const
  {
    x = feasible(x, data.M);
    y = feasible(y, data.N);

    return data(x,y);
  }

  int feasible(int x, int M) const
  {
    if (x < 0)
      return 0;

    if (x > M-1)
      return M-1;

    return x;
  }

protected:
	std::vector<R> x;
	std::vector<R> y;
	std::vector<int> linear_indices;

	const matrix<double> data;
	const vector<double> vd;
};

class Boundary_points_calculator
{
	public:
	  Boundary_points_calculator(const matrix<double>& data, const vector<double>& vd)
	    :  data(data), vd(vd), M(data.M), N(data.N)
	 {}

	 // Triplets
	 template<typename R>
	 Boundary_points<R> operator () (R x0, R y0, R x1, R y1, R x2, R y2) const 
	 {
	 		Boundary_points<R> points(data, vd, x0, y0);
	 		add_points_along_a_line_segment(x0, y0, x1, y1, points);
	 		add_points_along_a_line_segment(x1, y1, x2, y2, points);

	 		return points;
	 }

	 // Pair
	 template<typename R>
	 Boundary_points<R> operator () (R x0, R y0, R x1, R y1) const 
	 {
	 		Boundary_points<R> points(data, vd, x0, y0);
			add_points_along_a_line_segment(x0, y0, x1, y1, points);

			return points;
	 }

	// This is similar to evaluate_line_integral
	// For pair of points
	template<typename R>
	void add_points_along_a_line_segment(R sx, R sy, R ex, R ey, Boundary_points<R>& points) const
	{
		using std::abs;
		using std::sqrt;

		R dx = (ex - sx)*vd[0];
		R dy = (ey - sy)*vd[1];
		R line_length = sqrt(dx*dx + dy*dy);

		typedef std::pair<R, int> crossing;

		// Invariant of direction
		auto intersection =
		[&]
		(int index_change, R line_length, R dx, R start,
		 std::vector<crossing>& crossings) -> void
		{
			R k = dx/line_length;

			if (abs(k) < 1e-6)
				return;

			R current;
			if (k > 0)
			{
				current = 1.0 - fractional_part(start);
			}
			else
			{
				current = fractional_part(start);
				index_change = -index_change;
			}

			// -eps to avoid adding last point
			for (; current < abs(dx) - 1e-6; current = current +1)
				crossings.push_back( crossing( abs(current/dx) , index_change) );
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

		int source_id =  int( spii::to_double(sx) )  + int( spii::to_double(sy) )*M;

		scratch_space->clear();
		scratch_space->push_back( crossing(R(-1e-100), 0) );
		intersection(1,line_length, dx, sx, *scratch_space);
		intersection(M,line_length, dy, sy, *scratch_space);
		std::sort(scratch_space->begin(), scratch_space->end());

		// Remove small increments
		auto prev = scratch_space->begin();
		auto next = scratch_space->begin();
		next++;
		for (; next != scratch_space->end(); prev++, next++)
		{
			// (Sorted)
			if ( (next->first - prev->first) > 1e-6)
			{				
				R xc = sx + (next->first)*dx;
				R yc = sy + (next->first)*dy;
				points.add(xc,yc, source_id);
			}

			source_id += next->second;
		}

		points.add(ex,ey, source_id);
	}

  const matrix<double> data;
	const vector<double> vd;
	int M; // rows
	int N; // cols
};