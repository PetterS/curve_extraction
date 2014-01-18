% This script shown an example where the cost for each edge
% is explicitly given instead of interpolated from an image/volume.
addpath([fileparts(mfilename('fullpath')) filesep 'library']);

rng(0);
linear_data = rand(50,50);

% Define local connectivity
dims = ndims(linear_data);
radius = 1; % 4-conn
connectivity = get_all_directions(radius,dims);

% Simple edge cost: 1/2 current + 1/2 target voxel value
data = inf([size(linear_data) size(connectivity,1)]);
problem_size = size(linear_data);

for x = 1:size(data,1);
	for y = 1:size(data,2)
		for c= 1:size(connectivity,1);		
			tx = x + connectivity(c,1);
			ty = y + connectivity(c,2);
			
			if ((tx < 1) || (tx > problem_size(1)))
				continue;
			end
			
			if ((ty < 1) || (ty > problem_size(2)))
				continue;
			end
			
			data(x,y,c) = (linear_data(x,y)+ linear_data(tx,ty))/2;

		end
	end
end

% inf not implemented yet
data(data > 1e10) = 1e10;

%
start_set = false(problem_size);
end_set = false(problem_size);

start_set(:,1) = true;
end_set(:,end) = true;

%% Compare to linear interpolation
figure(1);
C = Curve_extraction('edge', data, connectivity, start_set, end_set);
C.solve();
C.display();
C1 = C.curve;

%% The connectivity is choosen such that the edge cost will be exactly the same
% as the explicitly given datacosts.
figure(2);
Cl = Curve_extraction('linear_interpolation', linear_data, start_set, end_set);
Cl.set_connectivity_by_radius(radius);
Cl.solve();
Cl.display();

C2 = C.curve;
assert(all(C1(:) == C2(:)));
