addpath([fileparts(mfilename('fullpath')) filesep 'library']);

rng(0);
linear_unary = rand(50,50);

% Define local connectivity
dims = ndims(linear_unary);
radius = 1; % 4-conn
connectivity = get_all_directions(radius,dims);

% Simple edge cost: 1/2 current + 1/2 target voxel value
unary = inf([size(linear_unary) size(connectivity,1)]);
problem_size = size(linear_unary);

for x = 1:size(unary,1);
	for y = 1:size(unary,2)
		for c= 1:size(connectivity,1);		
			tx = x + connectivity(c,1);
			ty = y + connectivity(c,2);
			
			if ((tx < 1) || (tx > problem_size(1)))
				continue;
			end
			
			if ((ty < 1) || (ty > problem_size(2)))
				continue;
			end
			
			unary(x,y,c) = (linear_unary(x,y)+ linear_unary(tx,ty))/2;

		end
	end
end

% inf not implemented yet
unary(unary > 1e10) = 1e10;

%
start_set = false(problem_size);
end_set = false(problem_size);

start_set(:,1) = true;
end_set(:,end) = true;

%% Compare to linear interpolation
figure(1);
C = Curve_extraction('edge', unary, connectivity, start_set, end_set);
C.solve();
C.display();

%%
figure(2);
Cl = Curve_extraction('linear_interpolation', linear_unary, start_set, end_set);
Cl.set_connectivity_by_radius(radius);
Cl.solve();
Cl.display();