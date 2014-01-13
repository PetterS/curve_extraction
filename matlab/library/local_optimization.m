% Wrapper for mex function.
function [optimized_path, solution_cost, solve_time] = local_optimization(mesh_map, data, input_path, settings)

problem_size = size(data);

% Check file modification dates and recompile mex file
my_name = mfilename('fullpath');
[base_path, base_name, ~] = fileparts(my_name);
addpath(base_path);

compile(base_path, base_name)

% Checks for disallowed pixels
data(mesh_map == 0) = inf;

 if length(problem_size) == 2
   input_path = [input_path ones(size(input_path,1),1)]; 
end

[optimized_path, solution_cost, solve_time] = local_optimization_mex(data, input_path, settings);

% Check that the solution fulfills constraints
assert( size(optimized_path,2) == size(input_path, 2) );

% Remove padding
if length(problem_size) == 2
	optimized_path	= optimized_path(:,1:2);
end
