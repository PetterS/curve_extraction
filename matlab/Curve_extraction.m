%% Class handling curve extraction
% For example usage, there are several examples in this folder.
%
% This class solves the shortest path problem.
% Data cost, start/end/disallowed sets are all restricted to be either 2D or 3D matrices of the same size.
% Right now there are three different types of problem formulation supported.
%
% data_type = 'linear_interpolation':
% The data cost of each edge is given by linearly interpolating the data matrix.
%
% data_type = 'edge'
% The data cost of each edge is given explicitly.
%
% data_type = 'geodesic'
% There is no data cost, the data matrix defines a height at each voxel.
% The shortest path is calculated on this surface.
%
classdef Curve_extraction < handle

	% Settings
	properties

		%%
		% The cost function for a give curve is defined as:
		% cost =  data_cost(curve)
		% 			+ length_penalty | length of curve |
		% 			+ curvature_penalty | curvature of curve |^curvature_power
		% 			+ torsion_penalty | torsion of curve|^torsion_power
		length_penalty = 0;
		curvature_penalty = 0,
		torsion_penalty = 0,
		curvature_power = 2.0;
		torsion_power = 2.0;

		%%
		% The above cost function penalizes either length, curvature or torsion.
		% The limits below are hard constraints on the curve.
		% Note: These are not upper  bounds on the cost function
		% but upper bounds on the length, curvature and torsion of the curve.
		% Modifying these results in a Constrained Shortest Path First (CSPF).
		length_limit = inf;
		curvature_limit = inf;
		torsion_limit = inf;

		% The CSPF is solved user supergradient optimization with the following settings:
		supergradient_max_duality_gap = 1e-8;
		supergradient_iterations = 25;

		% The range the supergradient optimizer searches over:
		% These settings should work for most applications.
		length_range = [0 1e6];
		curvature_range = [0 1e6];
		torsion_range = [0 1e6];

		% Override penalty functions to make corresponding solution
		% to obey he limits.
		% This is useful for local optimization post-processing.
		override_penalty = true;

		% When curvature regularization is used potential speedups can be achieved
		% by A*.
		use_a_star = true;

		% Verbose information from the optimization.
		verbose = false;

		% Store the visit order when performing shortest path calculations.
		store_visit_time = false;

		% Number of threads the optimizer uses.
		% Three parts of the code are parallelized.
		% 1. All costs in a neighborhood.
		% 2. Linearly interpolation a data cost.
		% 3. Local optimization.
		num_threads = int32(min(2,feature('numThreads')));

		% Data to optimize over
		data =  [];

		% Defines the connectivity as a delta functions
		% Each row is a new edge
		% E.g. [dx dy dz] for voxel (x,y,z) gives and edge to (x+dx,y+dy,z+dz)
		% Initial connectivity is set in the constructor.
		% Tip: Can be generated via the method obj.set_connectivity_by_radius.
		connectivity = [];

		% Settings local optimization
		local_optimization_maxiter = 1000;
		
		% When performing local optimization addtional points are interpolated
		% s.t. each curve segment is at most this long
		local_optimzation_max_curve_segment_length = 1;
		descent_method = 'lbfgs';

		% The length of each voxel side
		% [x_length y_length z_length]
		voxel_dimensions = [1 1 1];

		% Placeholder for a shortest path
		curve = [];

		% "Virtual" parameters set in the constructor.
		start_set;
		end_set;
		disallowed_set;
	end

	% Stored by the solver
	properties (SetAccess = protected)
		% In case of (CSPF)
		duality_gap = 0;

		data_type = '';
		cost;
		info;
		visit_map = [];
	end

	properties (Hidden)
		% Internal bookkeeping.
		problem_size;
		cached_cost;
		cached_info;
		mesh_map = [];
		performed_local_optimization = false;

		% Default settings when methods are not called with an explicit argument.
		default_connectivity_radius = 3;
	end

	methods (Access = protected)
		% Convert settings to format supported by the mex file.
		function settings = gather_settings(self)
			settings.length_penalty = self.length_penalty;
			settings.curvature_penalty = self.curvature_penalty;
			settings.torsion_penalty = self.torsion_penalty;
			settings.curvature_power = self.curvature_power;
			settings.use_a_star = self.use_a_star;
			settings.verbose = self.verbose;
			settings.data_type = self.data_type;
			settings.maxiter = self.local_optimization_maxiter;
			settings.store_visit_time =  self.store_visit_time;
			settings.descent_method = self.descent_method;
			settings.voxel_dimensions = self.voxel_dimensions;
			settings.num_threads = self.num_threads;
		end

		% Used when the shortest path on a surface is calculated.
		function geodesic_constructor(self, data, varargin)

			% Default length penalty is ~= 0 for this data type
			self.length_penalty = 1;

			if ~ismatrix(data)
				error('Only 2D dimensional data supported by geodesic shortest path');
			end

			self.interpolation_constructor(data, varargin{:});
		end

		% Used when the data cost is defined by by performing linear interpolation a 2D or 3D matrix.
		function interpolation_constructor(self, data, varargin)
			self.problem_size = size(data);
			self.data = data;

			self.set_connectivity_by_radius(self.default_connectivity_radius);

			self.create_mesh_map(varargin{:});
		end

		% Used when each edge in the graph is given an explicit cost.
		function edge_constructor(self, data, connectivity, varargin)
			self.connectivity = int32(connectivity);
			self.problem_size = size(data);
			self.data = data;

			sz = size(data);
			self.problem_size = sz(1:end-1);

			self.create_mesh_map(varargin{:});
		end

		% Create a mesh_map which is used internally to keep track on
		% start set, end_set and disallowed pixel/voxels.
		function self = create_mesh_map(self, start_set, end_set, disallowed)
			% Create a structure holding the start/end and allowed pixels
			% for the algorithm to visit.
			mesh_map = ones(self.problem_size,'uint8');

			if nargin > 3
				assert(all(self.problem_size == size(disallowed)));
				mesh_map(disallowed) = 0;
			end

			assert(all(self.problem_size == size(start_set)));
			mesh_map(start_set) = 2;

			% End set
			if nargin > 2
				if (any(start_set(:) & end_set(:)))
					error('Some voxels are both in the start and end set');
				end

				assert(all(self.problem_size == size(end_set)));
				mesh_map(end_set) = 3;
			end

			% Save
			self.mesh_map = mesh_map;
		end

		function feasible = curve_feasibility_check(self, info)
				length_check = info.length < self.length_limit;
				curvature_check = info.curvature < self.curvature_limit;
				torsion_check = info.torsion < self.torsion_limit;

				% Avoid potential numerical issues
				if (self.length_limit == inf)
					length_check = true;
				end
				if (self.curvature_limit == inf)
					curvature_check = true;
				end
				if (self.torsion_limit == inf)
					torsion_check = true;
				end

				feasible = length_check && curvature_check && torsion_check;
		end
	end
	methods

		% Constructor
		function self = Curve_extraction(data_type, varargin)
			addpath([fileparts(mfilename('fullpath')) filesep 'library']);

			self.data_type = data_type;

			% Send the input to the correct constructor:
			switch data_type
				case 'linear_interpolation'
					self.interpolation_constructor(varargin{:});
				case 'edge'
					self.edge_constructor(varargin{:});
				case 'geodesic'
					self.geodesic_constructor(varargin{:});
			end

			if (length(self.problem_size) ~= 2) && (length(self.problem_size) ~= 3)
				error('Only two- and three-dimensional problems supported.');
			end

		end

		% Shortest path from the start set to the end set.
		function  [curve, cost, time, evaluations, visit_map ] = shortest_path(self)
			if nargin < 2
				settings = gather_settings(self);
			end

			if (~any(self.mesh_map(:) == 3))
				error('The problem has no end set.');
			end

			[curve, ~, time, evaluations, visit_map] = ...
				curve_segmentation(self.data_type, self.mesh_map, self.data, self.connectivity, settings);

			% Saving solution
			self.curve = curve;
			self.visit_map = visit_map;
			self.performed_local_optimization = false;

			% When self.curve is set the cost is automatically updated with detailed info
			cost = self.cost;

			%% Supergradient
			% Check if curve is within limits otherwise run super-gradient optimization
			if (self.curve_feasibility_check(self.info))
				return;
			end

			% Variables
			% x(1): function value
			% x(2): length penalty variable
			% x(3): curvature penalty variable
			% x(4): torsion penalty variable
			c = [-1; 0; 0; 0];
			A = zeros(0, 4);
			b = zeros(0, 1);
			lb = [-inf; self.length_range(1); self.curvature_range(1); self.torsion_range(1)];
			ub = [inf;  self.length_range(2); self.curvature_range(2); self.torsion_range(2)];

			projected_cost = nan;

			% Have we found any feasible solution?
			feasible_solution = false;

			% Determine which constraints are needed
			% Constraint 1: Length
			% Constraint 2: Curvature
			% Constraint 3: Torsion
			active_constraints = true(3,1);
			if (self.length_limit == inf)
				active_constraints(1) = false;
			end

			if (self.curvature_limit == inf)
				active_constraints(2) = false;
			end

			if (self.torsion_limit == inf)
				active_constraints(3) = false;
			end

			initial_multipliers = [settings.length_penalty; settings.curvature_penalty; settings.torsion_penalty];
			multipliers = initial_multipliers;

			limits = [self.length_limit; self.curvature_limit; self.torsion_limit];
			limits(~active_constraints) = 0;

			for iter = 1:self.supergradient_iterations

				if iter > 1
					% Determine which point on the dual function to evaluate.
					options = optimset('Display','none');
					x = linprog(c, A, b, [], [], lb, ub,[], options);

					projected_cost = x(1);

					multipliers = x(2:4);
					multipliers(~active_constraints) = initial_multipliers(~active_constraints);

					settings.length_penalty     = multipliers(1);
					settings.curvature_penalty	= multipliers(2);
					settings.torsion_penalty	  = multipliers(3);

					[sg.curve, ~, sg.time, sg.evaluations, sg.visit_map] = ...
						curve_segmentation(self.data_type, self.mesh_map, self.data, self.connectivity, settings);

					time = time+sg.time;
					evaluations = evaluations + sg.evaluations;

					[cost,info] = curve_info(self.data_type, self.data, sg.curve, self.connectivity, settings);
				else

					% Shortest path has already been run with initial penalty settings.
					cost = self.cost;
					info = self.info;
				end

				values(1,1) = info.length;
				values(2,1) = info.curvature;
				values(3,1) = info.torsion;

				% Constraint g(x) <= 0.
				g = values - limits;
				g(~active_constraints) = 0;

				% Dual function.
				% The solver returns  \min f(x) + \sum \multiplier * value
				% the dual function is defined as d(x) =  \min f(x) + \sum \multiplier*(value -limit)
				dual_val = cost.total - sum(multipliers.*limits);

				% Duality gap
				gap = projected_cost - dual_val;

				A = [A; 1, -g(:)'];
				b = [b; dual_val - sum(multipliers.*g)];

				% Feasibility check
				if (self.curve_feasibility_check(info))
					feasible_solution = true;

					% Store feasible solution
					curve = sg.curve;
					visit_map = sg.visit_map;
					self.duality_gap = gap;
				end

				if (self.verbose)
					fprintf('Iteration: %d \n', iter);

					if (self.curve_feasibility_check(info))
						fprintf('Current solution is feasible. \n');
					else
						fprintf('Current solution is infeasible. \n');
					end

					fprintf('Duality gap: %g \n', gap);
					fprintf('Length: %g max: %g	penalty: %g \n', values(1), limits(1), multipliers(1));
					fprintf('Curvature: %g max: %g	penalty: %g. \n',  values(2), limits(2), multipliers(2));
					fprintf('Torsion: %g max: %g	penalty: %g. \n',  values(3), limits(3), multipliers(3));
					fprintf('---\n');
					self.plot_curve(curve);
					drawnow();
				end

				if gap < self.supergradient_max_duality_gap
					break;
				end
			end

			if (~feasible_solution)
				warning('Unable to find a feasible solution. Returning optimal solution for the orginal cost function.');
			else
				if (self.override_penalty)
					self.length_penalty = settings.length_penalty;
					self.curvature_penalty = settings.curvature_penalty;
					self.torsion_penalty = settings.torsion_penalty;
				end

				self.curve = curve;
				self.visit_map = visit_map;
			end
		end


		% Compute distance to every point from the start set
		function [distances, curve] = compute_all_distances(self)
			settings = gather_settings(self);
      settings.store_distances = true;
      settings.compute_all_distances = true;

			[curve, ~, ~, ~, ~, ~, distances] = ...
			 		 curve_segmentation(self.data_type, self.mesh_map, self.data, self.connectivity, settings);

			if 	min(distances(self.disallowed_set))	 < 1e38
				error('Disallowed set should have infinte cost.');
			end
			distances(self.disallowed_set) = inf;

		end

		% Compute visitation_tree
		% Returns a matrix (named tree) where each entry contains the linear index to
		% the parent voxel when every pixel/voxel is visited.
    function  [tree, curve] = compute_visit_tree(self)
			settings = gather_settings(self);
      settings.store_parents = true;
      settings.compute_all_distances = true;

			[curve, ~, ~, ~, ~, tree] = ...
			 		 curve_segmentation(self.data_type, self.mesh_map, self.data, self.connectivity, settings);

			% If no curve is store this
			if isempty(self.curve)
				self.curve = curve;
			end
		end

		% Perform local optimization on a discrete shortest_path solution.
		% The user can perform local optimization on any curve
		% by setting obj.curve and then calling this method.
		function [curve,cost,time] = local_optimization(self)

			if strcmp(self.data_type,'linear_interpolation')
			 %ok
			elseif strcmp(self.data_type,'geodesic')
				%ok
			else
				error('Local optimization is only supported  for data_type ={linear_interpolation, geodesic}');
			end

			% Add more points first time
			if ~self.performed_local_optimization
				curve = interpolate_more_points(self, self.curve, self.local_optimzation_max_curve_segment_length);
				self.performed_local_optimization = true;
			else
				curve = self.curve;
			end

			if isempty(self.curve)
				fprintf('No curve stored, running the shortest_path \n');
				self.shortest_path()
			end

			if (~any(self.mesh_map(:) == 3))
				error('The problem has no end set.');
			end

			settings = gather_settings(self);
			[curve, ~, time, success] = local_optimization(	self.data_type, curve, self.mesh_map, ...
																											self.data, self.connectivity, settings);

			if (~success)
				if (self.verbose)
					warning('Unable to find a better local optima. Keeping the old solution.');
				end

				curve = self.curve;
				cost = self.cost;
			else
				% Saving solution
				self.curve = curve;
				cost = self.cost;
			end
		end

		% Show information about the stored shortest path and settings.
		function display(self)
			details(self);
			self.plot_curve(self.curve);
		end

		function plot_curve(self, curve)
			if nargin < 2
				curve = self.curve;
			end

			clf; hold on;
			msgs = {};

			if (~isempty(curve))
				msgs{end+1} = sprintf('Cost; total: %g data: %g, length: %g curvature: %g, torsion %g.', ...
					self.cost.total, self.cost.data, self.cost.length, self.cost.curvature, self.cost.torsion);
			end

			msgs{end+1} = sprintf('Cost function; data + %g|length| + %g|curvature|^{%g} + %g|torsion|^{%g}.', ...
					self.length_penalty, self.curvature_penalty, self.curvature_power, self.torsion_penalty, self.torsion_power);

			if (~isempty(curve))
				msgs{end+1} = sprintf('Curve info: Length: %g |curvature|^{%g}: %g  |torsion|^{%g}: %g', ...
					self.info.length,  self.curvature_power, self.info.curvature,  self.torsion_power, self.info.torsion);
			end

			if (strcmp(self.data_type,'linear_interpolation') || strcmp(self.data_type,'geodesic'))

				if (length(self.problem_size) == 2)
					cost_im = self.data;
					cost_im(self.mesh_map ~= 1) = -1;
					imagesc(double(cost_im))
					colormap gray(256); axis equal; axis off;
					axis ij;
				end

			elseif (strcmp(self.data_type,'edge'))

				% Just display start,end,allowed set.
				if (length(self.problem_size) == 2)
					imagesc(3-self.mesh_map);
					colormap(gray(4)); axis equal; axis off;
					axis ij;
					title({msgs{:}});
				end
			end

			if (~isempty(curve))
				title({msgs{:}});

				cmap = jet(3);

				% Draw the stored solution.
				if (length(self.problem_size) == 3)
					plot3(curve(:,1), curve(:,2), curve(:,3),'-r', 'linewidth', 2);

					plot3(curve(1,1),curve(1,2),curve(1,3),'ko', 'MarkerFaceColor',cmap(2,:), 'MarkerSize',5);
					plot3(curve(end,1),curve(end,2),curve(end,3),'ko','MarkerFaceColor', cmap(3,:), 'MarkerSize',5);

					legend('Curve','Start','End','Location', 'EastOutside');
					view(3)
				else

					plot(curve(:,2),curve(:,1),'r-' , 'linewidth',2)
					plot(curve(1,2),curve(1,1),'ko','MarkerFaceColor', cmap(2,:), 'MarkerSize',5);
					plot(curve(end,2),curve(end,1),'ko','MarkerFaceColor', cmap(3,:), 'MarkerSize',5);

					legend('Curve','Start','End','Location', 'EastOutside');
				end
			else
				fprintf('No solution stored, please run obj.shortest_path() \n');
			end
		end
		
		% Take a curve and add points s.t. each line segment is at most
		% max_length long
		function new_curve = interpolate_more_points(self, curve,max_length)
						
			new_curve = curve(1,:);
			for iter = 2:size(self.curve,1)
				s = self.curve(iter-1,:);
				e = self.curve(iter,:);
				v = e-s;
				
				d = norm(v);
				
				% points between s and e (including s,e)
				np = ceil(d/max_length)+1;
				kl = linspace(0,1,np);
				kl = kl(2:end);
				
				for k = kl
					new_curve = [new_curve; s+v*k];
				end
			end
		end
		
		%% Set functions
		% This function updates to connectivity to consist of
		% every edge with a unique angle whose length is <= radius.
		function set_connectivity_by_radius(self, radius)

			if strcmp(self.data_type,'edge')
				error('Cannot modify connectivity for Curve_extraction object with data_type: edge.');
			end

			% Generate connectivity
			connectivity = get_all_directions(radius, length(self.problem_size));
			connectivity = int32(connectivity);

			if size(connectivity,2) == 2
				connectivity = [connectivity zeros(size(connectivity,1),1)];
			end

			self.connectivity = connectivity;
		end

		function set.num_threads(self, num_threads)
			assert(num_threads > 0);
			self.num_threads = int32(num_threads);
		end

		function set.voxel_dimensions(self, voxel_dimensions)
			assert(length(voxel_dimensions) == 2 || length(voxel_dimensions) == 3);
			assert(all(voxel_dimensions) > 0);

			if (length(voxel_dimensions) == 2)
				voxel_dimensions =[voxel_dimensions 1];
			end

			self.voxel_dimensions = voxel_dimensions;
			self.reset_solution();
		end

		% Two different local optimization algorithms are supported:
		% lbgfs (Broyden–Fletcher–Goldfarb–Shanno) uses symbolic differentiation.
		% nelder-mead: purely uses function evaluations and performs no differentiation.
		%
		% Tip: It is recommended to use lbgfs.
		function set.descent_method(self, method)
				switch method
					case 'lbfgs'
						%ok
					case 'nelder-mead'
						%ok
					otherwise
						error('Possible descent methods = {lbfgs, nelder-mead}');
				end

				self.descent_method = method;
		end

		function set.connectivity(self, connectivity)

			assert((size(connectivity,2) ~= 3) ||  (size(connectivity,2) ~= 2))

			% Padding
			if (size(connectivity,2) == 2)
				connectivity = [connectivity zeros(size(connectivity,1),1)];
			end

			% No duplicates
			if (numel(unique(connectivity,'rows')) ~= numel(connectivity))
				error('Duplicate entries in the given connectivity');
			end


			assert(size(connectivity,2) == 3);

			self.connectivity = connectivity;
		end

		function set.data(self, data)
			if (~isa(data,'double'));
				disp('Data-term must be a double, converting.');
				self.data = double(data);
			else
				self.data = data;
			end

			self.reset_solution();
		end

		function set.length_penalty(self, length_penalty)
			if (length_penalty < 0)
				error('Regularization coefficients must non-negative');
			end

			self.length_penalty = length_penalty;
			self.reset_solution();
		end

		function set.curvature_penalty(self, curvature_penalty)
			if (curvature_penalty < 0)
				error('Regularization coefficients must non-negative');
			end

			self.curvature_penalty = curvature_penalty;
			self.reset_solution();
		end

		function set.torsion_power(self, torsion_power)
			if (torsion_power  < 0)
				error('Torsion power must non-negative');
			end

			self.torsion_power = torsion_power;
			self.reset_solution();
		end

		function set.curvature_power(self, curvature_power)
			if (curvature_power  < 0)
				error('Torsion power must non-negative');
			end

			self.curvature_power = curvature_power;
			self.reset_solution();
		end

		function set.verbose(self, verbose)
			if ~isa(verbose,'logical')
				error('verbose must either be true or false');
			end

			self.verbose = verbose;
		end

		function set.store_visit_time(self, store_visit_time)
			if ~isa(store_visit_time,'logical')
				error('store_visit_time must either be true or false');
			end

			self.store_visit_time = store_visit_time;
		end

		function set.use_a_star(self, use_a_star)
			if ~isa(use_a_star,'logical')
				error('use_a_star must either be true or false');
			end

			self.use_a_star = use_a_star;
		end

		function set.torsion_penalty(self, torsion_penalty)
			if (torsion_penalty < 0)
				error('Regularization coefficients must non-negative');
			end

			self.torsion_penalty = torsion_penalty;
			self.reset_solution();
		end

		% Calculate curve cost on demand
		function cost = get.cost(self)
			if isempty(self.cached_cost)
				[self.cached_cost,self.cached_info] = self.curve_info(self.curve);
			end

			cost = self.cached_cost;
		end

		function info = get.info(self)
			if isempty(self.cached_cost)
				[self.cached_cost,self.cached_info] = self.curve_info(self.curve);
			end

			info = self.cached_info;
		end

		% Keeping the cuvre
		function reset_solution(self)
			self.cached_cost = [];
		end

		function set.data_type(self, data_type)

			switch data_type
				case 'linear_interpolation'
					%ok
				case 'edge'
					%ok
				case 'geodesic'
					%ok
				otherwise
					error(sprintf('Possible data types = {linear_interpolation, edge, geodesic}'));
			end

			self.data_type = data_type;
		end

		% Interface the start,end and disallowed set
		function allowed = get.disallowed_set(self)
			allowed = self.mesh_map == 0;
		end

		function set.disallowed_set(self, disallowed)
			assert(all(self.problem_size == size(disallowed)));

			self.mesh_map(self.mesh_map == 0) = 1;
			self.mesh_map(disallowed) = 0;
		end

		function start_set = get.start_set(self)
			start_set = self.mesh_map == 2;
		end

		function set.start_set(self, start_set)
			assert(all(self.problem_size == size(start_set)));

			remove = (self.mesh_map == 2) & ~start_set;
			self.mesh_map(start_set) = 2;
			self.mesh_map(remove) = 1;
		end

		function end_set = get.end_set(self)
			end_set = self.mesh_map == 3;
		end

		function set.end_set(self, end_set)
			assert(all(self.problem_size == size(end_set)));

			remove = (self.mesh_map == 3) & ~end_set;
			self.mesh_map(end_set) = 3;
			self.mesh_map(remove) = 1;
		end

		function set.mesh_map(self, mesh_map)

			if (~any(mesh_map(:) == 2))
				error('The problem has no start set.');
			end

			self.mesh_map = mesh_map;
			self.reset_solution();
		end

		function set.curve(self, curve)
			if (isempty(curve))
				return;
			end

			assert(size(curve,2) == length(self.problem_size));
			self.curve = curve;
			self.reset_solution();
		end

		function check_limit(~, limit)
			if (limit(1) < 0)
				error('Limit must be >= 0.');
			end
		end

		function set.length_limit(self, limit)
			self.check_limit(limit);

			self.length_limit = limit;
		end

		function set.curvature_limit(self, limit)
			self.check_limit(limit);

			self.curvature_limit = limit;
		end

		function set.torsion_limit(self, limit)
			self.check_limit(limit);

			if (numel(self.problem_size) == 2)
				error('2D problem always have zero torsion.');
			end
			self.torsion_limit = limit;
		end

		function check_range(~, range)
			if (numel(range) ~= 2)
				error('Range must be given as a vector as [lower_bound upper_bound].');
			end

			if (range(1) < 0)
				error('Lower bound must be >= 0.');
			end

			if (range(2) < range(1))
				error('upper bound < lower bound.');
			end
		end

		function set.length_range(self, range)
			self.check_range(range);

			self.length_range = range;
		end

		function set.curvature_range(self, range)
			self.check_range(range);

			self.curvature_range = range;
		end

		function set.torsion_range(self, range)
			self.check_range(range);

			self.torsion_range = range;
		end

		function set.local_optimization_maxiter(self, iterations)
			if (iterations < 0)
				error('Iterations >= 0.');
			end

			self.local_optimization_maxiter = iterations;
		end

		function set.supergradient_iterations(self, iterations)
			if (iterations < 0)
				error('Iterations >= 0.');
			end

			self.supergradient_iterations = iterations;
		end

		function set.supergradient_max_duality_gap(self, rel_gap)
			if (rel_gap < 0)
				error('Relative gap  >= 0.');
			end

			self.supergradient_max_duality_gap = rel_gap;
		end

		function set.local_optimzation_max_curve_segment_length(self, length)
			if (length <= 0)
				error('Length > 0')
			end
			
			self.local_optimzation_max_curve_segment_length = length;
		end
		
		% Given a curve this function calculates the total cost
		% of the curve given the current parameters and data cost
		% and stores it in cost.
		% Info contains the total, length, curvature and torsion
		% of the calculated with current curvature_power and torsion_power.
		function [cost,info] = curve_info(self, curve)

			if nargin < 2
				curve = self.curve;
			end

			if ~isempty(curve)
				settings = gather_settings(self);
				[cost,info] = curve_info(self.data_type, self.data, curve, self.connectivity, settings);
			else
				cost = [];
				info = [];
			end
		end
	end
end