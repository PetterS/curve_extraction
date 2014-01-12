// This function calculates the unary cost and 
//the cost of  length, curvature and torsion regularization. 
#include "curve_segmentation.h"

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ASSERT(nlhs == 5);
	ASSERT(nrhs == 3);

	// Parse data
	int curarg = 0;
	const matrix<double> unary_matrix(prhs[curarg++]);
	const matrix<double> path(prhs[curarg++]);
	MexParams params(nrhs-curarg, prhs+curarg);
  InstanceSettings settings = parse_settings(params); 


  // Returns length, curvature and unary cost
  matrix<double> unary_cost(1);
  matrix<double> length_cost(1);
  matrix<double> curvature_cost(1);
  matrix<double> torsion_cost(1);
  matrix<double> total_cost(1);

  total_cost(0) = 0;
  unary_cost(0) = 0;
  length_cost(0) = 0;
  curvature_cost(0) = 0;
  torsion_cost(0) = 0;

  plhs[0] = total_cost;
  plhs[1] = unary_cost;
  plhs[2] = length_cost;
  plhs[3] = curvature_cost;
  plhs[4] = torsion_cost;

  PieceWiseConstant data_term( unary_matrix.data,
                               unary_matrix.M,
                               unary_matrix.N,
                               unary_matrix.O,
                               settings.voxel_dimensions);

  // Unary and length
  for (int k = 0; k < path.M-1; k++)
  {
    unary_cost(0) += data_term.evaluate_line_integral
                          (path(k+0,0) -1, path(k+0,1) -1, path(k+0,2) -1,
                           path(k+1,0) -1, path(k+1,1) -1, path(k+1,2) -1);
  }

  if (settings.length_penalty > 0)
  {
    length_cost_functor length_cost_fun(settings.voxel_dimensions, settings.length_penalty);

    for (int k = 0; k < path.M-1; k++)
    {
      length_cost(0) += length_cost_fun(path(k+0,0) -1, path(k+0,1) -1, path(k+0,2) -1,
                                        path(k+1,0) -1, path(k+1,1) -1, path(k+1,2) -1);
    }
  }

 	// Curvature
  if (settings.curvature_penalty > 0)
  {
    curvature_cost_functor curvature_cost_fun(settings.voxel_dimensions, settings.curvature_penalty, settings.curvature_power);

   	for (int k = 0; k < path.M-2; k++)
   	{
       curvature_cost(0) += 
       curvature_cost_fun((path(k+0,0)-1),
             				      (path(k+0,1)-1),
             				      (path(k+0,2)-1),
             				      (path(k+1,0)-1),
             				      (path(k+1,1)-1),
             				      (path(k+1,2)-1),
              			      (path(k+2,0)-1),
              			      (path(k+2,1)-1),
              			      (path(k+2,2)-1));
   	}
  }

  if (settings.torsion_penalty > 0)
  {
    torsion_cost_functor torsion_cost_fun(settings.voxel_dimensions, settings.torsion_penalty, settings.torsion_power);

    // Torsion
    for (int k = 0; k < path.M-3; k++)
    {
      torsion_cost(0) += 
      torsion_cost_fun((path(k+0,0)-1),
                       (path(k+0,1)-1),
                       (path(k+0,2)-1),
                       (path(k+1,0)-1),
                       (path(k+1,1)-1),
                       (path(k+1,2)-1),
                       (path(k+2,0)-1),
                       (path(k+2,1)-1),
                       (path(k+2,2)-1),
                       (path(k+3,0)-1),
                       (path(k+3,1)-1),
                       (path(k+3,2)-1));
    }
  }

  total_cost(0) = unary_cost(0) + length_cost(0) + curvature_cost(0);
}