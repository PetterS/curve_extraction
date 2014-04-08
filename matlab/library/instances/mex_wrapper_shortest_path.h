#pragma once

template<typename Data_cost, typename Pair_cost, typename Triplet_cost, typename Quad_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction(int            nlhs,     /* number of expected outputs */
                 mxArray        *plhs[],  /* mxArray output pointer array */
                 int            nrhs,     /* number of inputs */
                 const mxArray  *prhs[]   /* mxArray input pointer array */)
{
 char problem_type[1024];
 if (mxGetString(prhs[0], problem_type, 1024))
   throw runtime_error("First argument must be a string.");

  if (!strcmp(problem_type,"linear_interpolation"))
    main_function< Linear_data_cost, Euclidean_length, Euclidean_curvature, Euclidean_torsion>(nlhs, plhs, nrhs, prhs);
  else if (!strcmp(problem_type,"edge"))
    main_function< Edge_data_cost, Euclidean_length, Euclidean_curvature, Euclidean_torsion>(nlhs, plhs, nrhs, prhs); 
  else if (!strcmp(problem_type,"geodesic"))
    main_function< Zero_data_cost, Geodesic_length, Zero_triplet, Zero_quad>(nlhs, plhs, nrhs, prhs);
  else
    throw runtime_error("Unknown data type");
}
