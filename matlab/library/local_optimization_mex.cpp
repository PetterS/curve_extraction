// Petter Strandmark and Johannes Ulén 2013.
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <limits>

#include <spii/auto_diff_term.h>
#include <spii/transformations.h>
#include <spii/solver.h>

// Mex
#include "curve_segmentation.h"

using namespace spii;
using namespace curve_extraction;

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

void mex_log_function(const std::string& str)
{
	mexPrintf("%s\n", str.c_str());
}

struct Point
{
	double xyz[3];
};

class Length
{
	public:
		Length(std::vector<double> dims_, double penalty_)
		{
			dims = dims_;
			penalty = penalty_;
		}

		template<typename R>
		R operator()(const R* const point1, const R* const point2) const
		{
			using std::sqrt;
			R dx = dims[0]*(point1[0] - point2[0]);
			R dy = dims[1]*(point1[1] - point2[1]);
			R dz = dims[2]*(point1[2] - point2[2]);

			return  penalty*sqrt(dx*dx + dy*dy + dz*dz);
		}

	private:
		std::vector<double> dims;
		double penalty;
};

class Curvature
{
public:
	Curvature(std::vector<double> dims_, double penalty_,  double power_)
	{
		dims = dims_;
		power = power_;
		penalty = penalty_;
	}

	template<typename R>
	R operator()(const R* const point1, const R* const point2, const R* const point3) const
	{
		return  penalty * compute_curvature<R>(
						point1[0]*dims[0], point1[1]*dims[1], point1[2]*dims[2],
						point2[0]*dims[0], point2[1]*dims[1], point2[2]*dims[2],
						point3[0]*dims[0], point3[1]*dims[1], point3[2]*dims[2],
						power);
	}

private:
	std::vector<double> dims;
	double power;
	double penalty;
};

class Torsion
{
public:
	Torsion(std::vector<double> dims_, double penalty_,  double power_)
	{
		dims = dims_;
		power = power_;
		penalty = penalty_;
	}

	template<typename R>
	R operator()(const R* const point1,
	             const R* const point2,
	             const R* const point3,
	             const R* const point4) const
	{
		return  penalty * compute_torsion<R>(
						point1[0]*dims[0], point1[1]*dims[1], point1[2]*dims[2],
						point2[0]*dims[0], point2[1]*dims[1], point2[2]*dims[2],
						point3[0]*dims[0], point3[1]*dims[1], point3[2]*dims[2],
						point4[0]*dims[0], point4[1]*dims[1], point4[2]*dims[2],
						power);
	}

private:
	std::vector<double> dims;
	double power;
	double penalty;
};

template<typename DataTermImpl>
class LinearUnary
{
public:
	LinearUnary(DataTermImpl& data_term_)
		: data_term(data_term_)
	{ }

	template<typename R>
	R operator()(const R* const point1, const R* const point2) const
	{
		return data_term.evaluate_line_integral(point1[0],point1[1],point1[2],
		                                        point2[0],point2[1],point2[2]);
	}

private:
	DataTermImpl& data_term;
};


void mexFunction_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double start_time = ::get_wtime();

	using namespace std;
	ASSERT(nrhs == 3);

	const matrix<double> data_matrix(prhs[0]);
	const matrix<double> path(prhs[1]);

	MexParams params(1, prhs+2);
	InstanceSettings settings = parse_settings(params);

	const int n = path.M;
	const int dim = path.N;

	ASSERT( (dim == 2) || (dim == 3) );

	if (settings.verbose)
	{
		mexPrintf("Solving using : %s \n", settings.descent_method_str.c_str());
		mexPrintf("Maximum iterations: %d \n", settings.maxiter);
		mexPrintf("function_improvement_tolerance: %g \n", settings.function_improvement_tolerance);
		mexPrintf("argument_improvement_tolerance: %g \n", settings.argument_improvement_tolerance);
	}

	// Function to be optimized
	Function f;

	// Create the points and add them as variables
	// to the function.
	vector<Point> points(n);

	for (int i = 0; i < n; i++)
	{
		// Zero based index
		points[i].xyz[0] = path(i,0) - 1;
		points[i].xyz[1] = path(i,1) - 1;
		points[i].xyz[2] = path(i,2) - 1;

		for (int j = 0; j < 3; ++j) {
			points[i].xyz[j] = max(points[i].xyz[j], 0.0);
		}


		points[i].xyz[0] = min(points[i].xyz[0], data_matrix.M - 1.0);
		points[i].xyz[1] = min(points[i].xyz[1], data_matrix.N - 1.0);
		points[i].xyz[2] = min(points[i].xyz[2], data_matrix.O - 1.0);

		f.add_variable(points[i].xyz, 3);
	}

	// The start and end of the curve is fixed.
	f.set_constant(points[0].xyz, true);
	f.set_constant(points[1].xyz, true);

	f.set_constant(points[n-2].xyz, true);
	f.set_constant(points[n-1].xyz, true);

	PieceWiseConstant data_term( data_matrix.data,
	                             data_matrix.M,
	                             data_matrix.N,
	                             data_matrix.O,
	                             settings.voxel_dimensions);

	// Adding data cost
	if (settings.data_type == linear_interpolation)
	{
		auto data = std::make_shared<AutoDiffTerm<LinearUnary<PieceWiseConstant>, 3, 3>>(data_term);

		for (int i = 1; i < n; ++i)
			f.add_term(data, points[i-1].xyz, points[i].xyz);
	}

	// Functor for each type of regularization penalty
	auto length = std::make_shared<AutoDiffTerm<Length, 3, 3>>(settings.voxel_dimensions, settings.length_penalty);
	auto curvature = std::make_shared<AutoDiffTerm<Curvature, 3, 3, 3>>
				(settings.voxel_dimensions, settings.curvature_penalty, settings.curvature_power);
	auto torsion = std::make_shared<AutoDiffTerm<Torsion, 3, 3, 3, 3>>
				(settings.voxel_dimensions, settings.torsion_penalty, settings.torsion_power);

	// Adding length cost
	if (settings.length_penalty > 0)
	{
		for (int i = 1; i < n; ++i)
			f.add_term(length, points[i-1].xyz, points[i].xyz);
	}

	// Add curvature cost
	if (settings.curvature_penalty > 0)
	{
		for (int i = 2; i < n; ++i)
		{
			vector<double*> args;
			args.push_back(points[i-2].xyz);
			args.push_back(points[i-1].xyz);
			args.push_back(points[i].xyz);

			f.add_term(curvature, args);
		}
	}

	// Adding torsion cost
	if (settings.torsion_penalty > 0)
	{
		for (int i = 3; i < n; ++i)
		{
			vector<double*> args;
			args.push_back(points[i-3].xyz);
			args.push_back(points[i-2].xyz);
			args.push_back(points[i-1].xyz);
			args.push_back(points[i].xyz);

			f.add_term(torsion, args);
		}
	}
	if (settings.verbose)
		mexPrintf("Initial function value: %.3e\n", f.evaluate());

	std::unique_ptr<Solver> solver;

  if (settings.descent_method == lbfgs) 
		solver.reset(new LBFGSSolver);
	else if (settings.descent_method == nelder_mead)
		solver.reset(new NelderMeadSolver);

	SolverResults results;

	if (settings.verbose)
		solver->log_function = mex_log_function;
	else
		solver->log_function = nullptr;
	

	solver->maximum_iterations = settings.maxiter;
	solver->function_improvement_tolerance = settings.function_improvement_tolerance;
	solver->argument_improvement_tolerance = settings.argument_improvement_tolerance;

	if (settings.num_threads > 0)
		f.set_number_of_threads(settings.num_threads);

	solver->solve(f, &results);

	if (settings.verbose)
	{
		std::stringstream sout;
		sout << results << endl;
		f.print_timing_information(sout);
		mexPrintf("%s\n", sout.str().c_str());
		mexPrintf("Final function value:   %.3e\n", f.evaluate());
	}


	matrix<double> resulting_path(points.size(), dim);
	for (int i = 0; i <n; i++)
	{
		resulting_path(i,0) = points[i].xyz[0] + 1;
		resulting_path(i,1) = points[i].xyz[1] + 1;
		resulting_path(i,2) = points[i].xyz[2] + 1;
	}

	plhs[0] = resulting_path;

	// Info from solver
	matrix<double> o_time(1);
	matrix<double> o_cost(1);

	double end_time = ::get_wtime();

	o_cost(0) = f.evaluate();
	o_time(0) = end_time - start_time;

	// Solution info
	plhs[1] = o_cost;
	plhs[2] = o_time;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	try {
		mexFunction_main(nlhs, plhs, nrhs, prhs);
	}
	catch (std::exception& error) {
		mexErrMsgTxt(error.what());
	}
	catch (...) {
		mexErrMsgTxt("Unknown exception.");
	}
}
