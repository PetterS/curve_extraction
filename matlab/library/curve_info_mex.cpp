// This function calculates the data cost and 
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
	ASSERT(nrhs == 4 || nrhs == 3);

	// Parse data
	int curarg = 0;
	const matrix<double> data_matrix(prhs[curarg++]);
	const matrix<double> path(prhs[curarg++]);
  const matrix<int> connectivity(prhs[curarg++]);
	MexParams params(nrhs-curarg, prhs+curarg);
  InstanceSettings settings = parse_settings(params); 


  // Returns length, curvature and data cost
  matrix<double> total_cost(1);

  matrix<double> total_length_cost(1);
  matrix<double> total_data_cost(1);
  matrix<double> total_curvature_cost(1);
  matrix<double> total_torsion_cost(1);
  

  total_cost(0) = 0;
  total_data_cost(0) = 0;
  total_length_cost(0) = 0;
  total_curvature_cost(0) = 0;
  total_torsion_cost(0) = 0;

  plhs[0] = total_cost;
  plhs[1] = total_data_cost;
  plhs[2] = total_length_cost;
  plhs[3] = total_curvature_cost;
  plhs[4] = total_torsion_cost;

  Data_cost data_cost(data_matrix, connectivity, settings);

  // Unary and length
  for (int k = 0; k < path.M-1; k++)
  {
    total_data_cost(0) += data_cost( path(k+0,0) -1, path(k+0,1) -1, path(k+0,2) -1,
                                path(k+1,0) -1, path(k+1,1) -1, path(k+1,2) -1);
  }

  if (settings.length_penalty > 0)
  {
    Length_cost length_cost_fun(settings.voxel_dimensions, settings.length_penalty);

    for (int k = 0; k < path.M-1; k++)
    {
      total_length_cost(0) += length_cost_fun(path(k+0,0) -1, path(k+0,1) -1, path(k+0,2) -1,
                                        path(k+1,0) -1, path(k+1,1) -1, path(k+1,2) -1);
    }
  }

 	// Curvature
  if (settings.curvature_penalty > 0)
  {
    Curvature_cost curvature_cost_fun(settings.voxel_dimensions, settings.curvature_penalty, settings.curvature_power);

   	for (int k = 0; k < path.M-2; k++)
   	{
       total_curvature_cost(0) += 
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
    Torsion_cost torsion_cost_fun(settings.voxel_dimensions, settings.torsion_penalty, settings.torsion_power);

    // Torsion
    for (int k = 0; k < path.M-3; k++)
    {
      total_torsion_cost(0) += 
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

  total_cost(0) = total_data_cost(0) + total_length_cost(0) + total_curvature_cost(0) + total_torsion_cost(0);
}