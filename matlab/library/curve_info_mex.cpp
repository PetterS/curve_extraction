// This function both calculates the cost of a curve decomposed into
// weighted and non weighted data, pair cost, triplet cost, quadruple cost
#include "curve_segmentation.h"

// Calls main_function
#include "instances/mex_wrapper_shortest_path.h"

template<typename Data_cost, typename Pair_cost, typename Triplet_cost, typename Quad_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ASSERT(nlhs == 8);
	ASSERT(nrhs == 4 || nrhs == 5);

	// Parse data
	int curarg = 1;
	const matrix<double> data_matrix(prhs[curarg++]);
	const matrix<double> input_path(prhs[curarg++]);
  const matrix<int> connectivity(prhs[curarg++]);
	MexParams params(nrhs-curarg, prhs+curarg);
  InstanceSettings settings = parse_settings(params);
  InstanceSettings info_settings = settings;


  // Pad to 3Ds
  matrix<double> path(input_path.M,3);
  for (int i = 0; i < input_path.M; i++)
  {
     path(i,0) = input_path(i,0);
     path(i,1) = input_path(i,1);

     if (input_path.N == 2)
      path(i,2) == 1;
     else
      path(i,2) = input_path(i,2); 
  }

  for (double& penalty : info_settings.penalty)
    penalty = 1;

  matrix<double> total_cost(1);

  matrix<double> total_data_cost(1);
  matrix<double> total_pair_cost(1);
  matrix<double> total_triplet_cost(1);
  matrix<double> total_quadruplet_cost(1);
  
  matrix<double> curve_pair(1);
  matrix<double> curve_triplet(1);
  matrix<double> curve_quadruplet(1);

  // Weighted by penalty function
  total_cost(0) = 0;
  total_data_cost(0) = 0;
  total_pair_cost(0) = 0;
  total_triplet_cost(0) = 0;
  total_quadruplet_cost(0) = 0;

  // E.g. Length, curvature and torsion of curve
  curve_pair(0) = 0;
  curve_triplet(0) = 0;
  curve_quadruplet(0) = 0;

  plhs[0] = total_cost;
  plhs[1] = total_data_cost;
  plhs[2] = total_pair_cost;
  plhs[3] = total_triplet_cost;
  plhs[4] = total_quadruplet_cost;
  plhs[5] = curve_pair;
  plhs[6] = curve_triplet;
  plhs[7] = curve_quadruplet;

  Data_cost data_cost(data_matrix, connectivity, info_settings);
  Pair_cost pair_cost(data_matrix, info_settings);

  // E.g. Data cost and length cost
  for (int k = 0; k < (int)path.M-1; k++)
  {
    Point p1(path(k+0,0) -1, path(k+0,1) -1, path(k+0,2) -1);
    Point p2(path(k+1,0) -1, path(k+1,1) -1, path(k+1,2) -1);

    total_data_cost(0)  += data_cost(  p1.xyz, p2.xyz);
    curve_pair(0)       += pair_cost(  p1.xyz, p2.xyz);
  }

 	// E.g curvature
  Triplet_cost triplet_cost(data_matrix, info_settings);
 	for (int k = 0; k < (int)path.M-2; k++)
 	{
    Point p1(path(k+0,0) -1, path(k+0,1) -1, path(k+0,2) -1);
    Point p2(path(k+1,0) -1, path(k+1,1) -1, path(k+1,2) -1);
    Point p3(path(k+2,0) -1, path(k+2,1) -1, path(k+2,2) -1);
   
    curve_triplet(0) += triplet_cost(p1.xyz, p2.xyz, p3.xyz);
 	}

  // E.g. Torsion
  Quad_cost quad_cost(data_matrix, info_settings);
  for (int k = 0; k < (int)path.M-3; k++)
  {
    Point p1(path(k+0,0) -1, path(k+0,1) -1, path(k+0,2) -1);
    Point p2(path(k+1,0) -1, path(k+1,1) -1, path(k+1,2) -1);
    Point p3(path(k+2,0) -1, path(k+2,1) -1, path(k+2,2) -1);
    Point p4(path(k+3,0) -1, path(k+3,1) -1, path(k+3,2) -1);

    curve_quadruplet(0) += quad_cost(p1.xyz, p2.xyz, p3.xyz, p4.xyz);
  }


  // Add weights
  // The if statement is added in order to avoid potential
  // numerical issues in the curvature/torsion calculations.
  if (settings.penalty[0] > 0)
    total_pair_cost(0) = curve_pair(0)*settings.penalty[0];
  else
    total_pair_cost(0) = 0;

  if (settings.penalty[1] > 0)
    total_triplet_cost(0) = curve_triplet(0)*settings.penalty[1];
  else
    total_triplet_cost(0) = 0;

  if (settings.penalty[2] > 0)
    total_quadruplet_cost(0) = curve_quadruplet(0)*settings.penalty[2];
  else
    total_quadruplet_cost(0) = 0;

  total_cost(0) = total_data_cost(0) + total_pair_cost(0) + total_triplet_cost(0) + total_quadruplet_cost(0);
}