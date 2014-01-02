% Wrapper for mex function.
function [curve, time, evaluations, cost, connectivity, visit_map ] ...
 = curve_segmentation(mesh_map, unary, settings, dirs)

% Parse option struct
if nargin >= 3
   
    % Only support struct input for simplicity
    assert(isa(settings,'struct'))

    % Converting to 0-based coordinate system used in C++ code
    % from 1-based coordinate system in matlab.
    if isfield(settings,'start_sets')
        for i = 1:length(settings.start_sets)
            settings.start_sets{i} = settings.start_sets{i} - 1;
        end
    end
    if isfield(settings,'end_sets')
        for i = 1:length(settings.end_sets)
            settings.end_sets{i} = settings.end_sets{i} - 1;
        end
    end
    
    % Default
    if ~isfield(settings,'voxeldimensions')
        settings.voxeldimensions = [1 1 1];
    end
    
    if ~ isfield(settings,'regularization_radius')
        settings.regularization_radius = 4.0;
    end
end

unary(mesh_map == 0) = max( max(unary(:))*1e2, 1e10);
assert(isa(mesh_map,'int32'));
assert(isa(unary,'double'));

if (~any(mesh_map(:) == 1))
    error('Atleast one pixel needs to be a normal visitable pixel.');
end

if ~any(mesh_map(:) == 2) && ~isfield(settings, 'start_sets')
    error('No start set');
end

if (~any(mesh_map(:) == 3)) && ~isfield(settings, 'end_sets')
    error('No end set');
end

assert(min(unary(:)) >= 0);
assert( all( size(mesh_map) == size(unary)));

% Check file modification dates and recompile mex file
my_name = mfilename('fullpath');
[base_path, base_name, ~] = fileparts(my_name);

extra_args{1} = ['-I' base_path];

sources{1} = 'node_segmentation.cpp'; % Length
sources{2} = 'edge_segmentation.cpp'; % Curvature
sources{3} = 'edgepair_segmentaion.cpp'; % Torsion

compile(base_path, base_name, sources, extra_args)


% If no connectivity is given one is generated 
% using regularization radius
if nargin < 4
	% Generate connectivity
	[~, dirs] = get_all_directions(settings.regularization_radius, ndims(unary));
	dirs = int32(dirs);
	
	if size(dirs,2) == 2
		dirs = [dirs zeros(size(dirs,1),1)];
	end	
end

[curve, time, evaluations, cost, connectivity, visit_map ] ...
 = curve_segmentation_mex(mesh_map, unary, dirs, settings);

%% Post process
% Remove third dim for 2 images
if (ndims(unary) == 2)
    curve = curve(:,1:2);
end

%return to coordinate system starting with 1.
curve = curve + 1;