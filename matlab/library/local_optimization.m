% Wrapper for mex function.
function [optimized_path, solution_cost, solve_time, success] =  ...
	local_optimization(	problem_type, input_path, mesh_map, ...
									 		data, connectivity, settings)

settings = parse_settings(settings);
problem_size = size(data);

% Check file modification dates and recompile mex file
my_name = mfilename('fullpath');
[base_path, base_name, ~] = fileparts(my_name);
addpath(base_path);

compile(base_path, base_name)

% Checks for disallowed pixels
data(mesh_map == 0) = inf;

if (~isa(data,'double'));
	disp('Data-term must be a double, converting.');
	data = double(data);
end

 if length(problem_size) == 2
   input_path = [input_path ones(size(input_path,1),1)]; 
end

problem_type = 'linear_interpolation';
[optimized_path, solution_cost, solve_time, success] ...
 = local_optimization_mex(problem_type, data, input_path, connectivity, settings);

% Check that the solution fulfills constraints
assert( size(optimized_path,2) == size(input_path, 2) );

% Remove padding
if length(problem_size) == 2
	optimized_path	= optimized_path(:,1:2);
end
