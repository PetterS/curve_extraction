% Wrapper for mex function.
function [curve, cost, time, evaluations, visit_map, shortest_path_tree] ...
 = curve_segmentation(problem_type, mesh_map, data, dirs, settings)

settings = parse_settings(settings);

if (~isa(mesh_map, 'uint8'))
    warning('mesh_map is not unsigned char (uint8) convering.')
    mesh_map = uint8(mesh_map);
end

data(mesh_map == 0) = inf;

if (~isa(mesh_map, 'uint8'))
    warning('mesh_map is not unsigned char (uint8) convering.')
    mesh_map = uint8(mesh_map);
end

if (~isa(data,'double'));
	disp('Data-term must be a double, converting.');
	data = double(data);
end

if (~any(mesh_map(:) == 1))
    error('Atleast one pixel needs to be a normal visitable pixel.');
end

if ~any(mesh_map(:) == 2) && ~isfield(settings, 'start_sets')
    error('No start set');
end

if (~any(mesh_map(:) == 3)) && ~isfield(settings, 'end_sets')
    error('No end set');
end

assert(min(data(:)) >= 0);

% Check file modification dates and recompile mex file
my_name = mfilename('fullpath');
[base_path, base_name, ~] = fileparts(my_name);
extra_args{1} = ['-I' base_path];
sources = {};

compile(base_path, base_name, sources, extra_args)

[curve, cost, time, evaluations, visit_map, shortest_path_tree] ...
 = curve_segmentation_mex(problem_type,  mesh_map, data, dirs, settings);

%% Post process
% Remove third dim for 2 images
if (ndims(mesh_map) == 2)
    curve = curve(:,1:2);
end

%return to coordinate system starting with 1.
curve = curve + 1;
shortest_path_tree = shortest_path_tree + 1;
