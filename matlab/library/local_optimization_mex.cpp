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

#include "curve_segmentation.h"

using namespace spii;
using namespace curve_extraction;

void mex_log_function(const std::string& str)
{
	mexPrintf("%s\n", str.c_str());
}

// Calls main_function
#include "instances/mex_wrapper_local_optimization.h"

template<typename Data_cost, typename Pair_cost, typename Triplet_cost, typename Quad_cost, typename Pentuple_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double start_time = ::get_wtime();

	using namespace std;
	ASSERT(nrhs == 5);
	ASSERT(nlhs == 4)

	int curarg = 1;
	const matrix<double> data_matrix(prhs[curarg++]);
	const matrix<double> input_path(prhs[curarg++]);
	const matrix<int> connectivity(prhs[curarg++]);

	MexParams params(nrhs-curarg, prhs+curarg);
	InstanceSettings settings = parse_settings(params);

	const int n = input_path.M;
	const int dim = input_path.N;

	ASSERT( (dim == 2) || (dim == 3) );

  // Pad to 3D
  matrix<double> path(n,3);
  for (int i = 0; i < n; i++)
  {
     path(i,0) = input_path(i,0);
     path(i,1) = input_path(i,1);

     if (dim == 2)
      path(i,2) == 1;
     else
      path(i,2) = input_path(i,2); 
  }


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
	 			(data_matrix, connectivity, settings);

	for (int i = 1; i < n; ++i)
		f.add_term(data, points[i-1].xyz, points[i].xyz);


	// Functor for each type of regularization penalty
	auto length = std::make_shared<AutoDiffTerm<Pair_cost, 3, 3>>
				(data_matrix, settings);
	auto curvature = std::make_shared<AutoDiffTerm<Triplet_cost, 3, 3, 3>>
				(data_matrix, settings);
	auto torsion = std::make_shared<AutoDiffTerm<Quad_cost, 3, 3, 3, 3>>
				(data_matrix, settings);
	auto jounce = std::make_shared<AutoDiffTerm<Pentuple_cost, 3, 3, 3, 3, 3>>
				(data_matrix, settings);


	if (settings.penalty[0] > 0)
	{
		for (int i = 1; i < n; ++i)
			f.add_term(length, points[i-1].xyz, points[i].xyz);
	}

	if (settings.penalty[1] > 0)
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

	if (settings.penalty[2] > 0)
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

	if (settings.penalty[3] > 0)
	{
		for (int i = 4; i < n; ++i)
		{
			vector<double*> args;
			args.push_back(points[i-4].xyz);
			args.push_back(points[i-3].xyz);
			args.push_back(points[i-2].xyz);
			args.push_back(points[i-1].xyz);
			args.push_back(points[i].xyz);

			f.add_term(jounce, args);
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
		for (int d = 0; d < dim; d++)
			resulting_path(i,d) = points[i].xyz[d] + 1;

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