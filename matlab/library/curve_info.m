% Wrapper for mex function.
function cost = curve_info(data, path, connectivity, settings)

% Check file modification dates and recompile mex file
my_name = mfilename('fullpath');
[base_path, base_name, ~] = fileparts(my_name);
addpath([base_path filesep '..']);

compile(base_path, base_name)

if (size(path,2) == 2)
   path = [path ones(size(path,1),1)]; 
end

[cost.total, cost.data, cost.length, cost.curvature, cost.torsion] ...
	= curve_info_mex(data, path, connectivity, settings);
