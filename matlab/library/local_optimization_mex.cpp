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

void mex_log_function(const std::string& str)
{
	mexPrintf("%s\n", str.c_str());
}

template<typename Data_cost, typename Length_cost, typename Curvature_cost, typename Torsion_cost>
void local_optimization(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double start_time = ::get_wtime();

	using namespace std;
	ASSERT(nrhs == 5);
	ASSERT(nlhs == 4)

	int curarg = 1;
	const matrix<double> data_matrix(prhs[curarg++]);
	const matrix<double> path(prhs[curarg++]);
	const matrix<int> connectivity(prhs[curarg++]);

	MexParams params(nrhs-curarg, prhs+curarg);
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

	// Adding data cost
	 auto data = std::make_shared<AutoDiffTerm<Data_cost, 3, 3>>
	 			(data_matrix, connectivity, settings.voxel_dimensions);

	for (int i = 1; i < n; ++i)
		f.add_term(data, points[i-1].xyz, points[i].xyz);


	// Functor for each type of regularization penalty
	auto length = std::make_shared<AutoDiffTerm<Length_cost, 3, 3>>
				(data_matrix, settings.voxel_dimensions, settings.length_penalty);
	auto curvature = std::make_shared<AutoDiffTerm<Curvature_cost, 3, 3, 3>>
				(data_matrix, settings.voxel_dimensions, settings.curvature_penalty, settings.curvature_power);
	auto torsion = std::make_shared<AutoDiffTerm<Torsion_cost, 3, 3, 3, 3>>
				(data_matrix, settings.voxel_dimensions, settings.torsion_penalty, settings.torsion_power);

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

	double intial_function_value  = f.evaluate();

	solver->maximum_iterations = settings.maxiter;
	solver->function_improvement_tolerance = settings.function_improvement_tolerance;
	solver->argument_improvement_tolerance = settings.argument_improvement_tolerance;

	if (settings.num_threads > 0)
		f.set_number_of_threads(settings.num_threads);

	solver->solve(f, &results);

	double optimized_function_value = f.evaluate();

	bool decreased_function_value = (optimized_function_value < intial_function_value)?true:false;

	if (settings.verbose)
	{
		std::stringstream sout;
		sout << results << endl;
		f.print_timing_information(sout);
		mexPrintf("%s\n", sout.str().c_str());

		if (decreased_function_value)
			mexPrintf("Local optimization successfully decreased the function value. \n");
		else
			mexPrintf("Local optimization unable to decrease the function value. \n");

		
		mexPrintf("Initial function value: %e\n", intial_function_value);
		mexPrintf("Final function value:   %e\n", optimized_function_value);	
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
	matrix<bool> o_successful(1);

	double end_time = ::get_wtime();

	o_cost(0) = optimized_function_value;
	o_time(0) = end_time - start_time;
	o_successful(0) = decreased_function_value;

	// Solution info
	plhs[1] = o_cost;
	plhs[2] = o_time;
	plhs[3] = o_successful;
}

// Wrapper data from MATLAB.
void mexFunction(int            nlhs,     /* number of expected outputs */
                 mxArray        *plhs[],  /* mxArray output pointer array */
                 int            nrhs,     /* number of inputs */
                 const mxArray  *prhs[]   /* mxArray input pointer array */)
{
 char problem_type[1024];
 if (mxGetString(prhs[0], problem_type, 1024))
   throw runtime_error("First argument must be a string.");

  if (!strcmp(problem_type,"linear_interpolation"))
    local_optimization<Linear_data_cost, Euclidean_length, Euclidean_curvature, Euclidean_torsion>(nlhs, plhs, nrhs, prhs);
  else if (!strcmp(problem_type,"geodesic"))
    local_optimization< Zero_data_cost, Geodesic_length, Geodesic_curvature, Zero_torsion>(nlhs, plhs, nrhs, prhs);
  else
    throw runtime_error("Unknown data type");
}