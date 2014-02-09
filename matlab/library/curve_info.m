% Wrapper for mex function.
function [cost, info] = curve_info(data, path, connectivity, settings)

settings = parse_settings(settings);

if (isempty(path))
	cost.total = nan;
	cost.curvature = nan;
	cost.data = nan;
	cost.length = nan;
	cost.torsion = nan;
	
	return;
end

if (~isa(data,'double'));
	disp('Data-term must be a double, converting.');
	data = double(data);
end

% Check file modification dates and recompile mex file
my_name = mfilename('fullpath');
[base_path, base_name, ~] = fileparts(my_name);
addpath([base_path filesep '..']);

compile(base_path, base_name)

if (size(path,2) == 2)
   path = [path ones(size(path,1),1)];
end

[cost.total, cost.data, cost.length, cost.curvature, cost.torsion, ...
 info.length, info.curvature, info.torsion] ...
	= curve_info_mex(data, path, connectivity, settings);
info.data = cost.data;