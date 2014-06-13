%% Class handling curve extraction
% For example usage see /examples/.
% This is the base class for all inherited versions.
%
% This class solves the shortest path problem.
% Data cost, start/end/disallowed sets are all restricted to be either 2D or 3D matrices of the same size.
%
%
classdef (Abstract) Curve_extraction_base < handle

	% To be set by a concrete subclass
	properties (Abstract, SetAccess = protected)
		cost;
		info;
		data_type;
	end

	properties
		% The CSPF is solved user supergradient optimization with the following settings:
		supergradient_max_duality_gap = 1e-8;
		supergradient_iterations = 25;

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
		num_threads = int32(1);

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
		% [x_length y_length z_length t_length]
		voxel_dimensions = ones(4,1);

		% Placeholder for a shortest path
		curve = [];

		% "Virtual" parameters set in the constructor.
		start_set;
		end_set;
		disallowed_set;
	end

	% In the inheriting classes make set and get methods
	% with appropriate name for the task.
	properties (Hidden)
		%%
		% The cost function for a give curve is defined as:
		% cost =  data_cost(curve)
		% 			+ penalty(n) | function(n) |^power(n)
		% where
		% function(n) = \inf if  function(n) > limit(n)
		%
		% Note: there is no "1 point penalty" the variable.
		% settings penalty(1) will not effect the results.
		%  
		penalty = zeros(4,1);
		power = ones(4,1);

		% These are local constraints on each quantity
		local_limit = inf(4,1);

		% These are global constraints on each quantity
		% The are hard constraints on the integral of the curve.
		% Note: These are not upper  bounds on the cost function
		% but upper bounds on the length, curvature and torsion of the curve.
		% Modifying these results in a Constrained Shortest Path (CSP).
		global_limit = inf(4,1)

		% The range the supergradient optimizer searches over:
		% These settings should work for most applications.
		range = [0 1e6;
						 0 1e6;
						 0 1e6];

		% Internal bookkeeping.
		problem_size;
		cached_cost;
		cached_info;
		mesh_map = [];
		checked_max_curve_segment_length = false;

		% Default settings when methods are not called with an explicit argument.
		default_connectivity_radius = 3;
	end

	% Stored by the solver
	properties (SetAccess = protected)
		% In case of (CSP)
		duality_gap = 0;
		visit_map = [];
	end


	methods (Access = protected)
		% Convert settings to format supported by the mex file.
		function settings = gather_settings(self)
			settings.penalty = self.penalty;
			settings.power = self.power;
			settings.local_limit = self.local_limit;

			settings.use_a_star = self.use_a_star;
			settings.verbose = self.verbose;
			settings.data_type = self.data_type;
			settings.maxiter = self.local_optimization_maxiter;
			settings.store_visit_time =  self.store_visit_time;
			settings.descent_method = self.descent_method;
			settings.voxel_dimensions = self.voxel_dimensions;
			settings.num_threads = self.num_threads;
			
			settings = self.parse_settings(settings);
		end
		
		% Make sure all data is in correct form.
		function settings = parse_settings(~, settings)
			int_entries = {'num_threads', 'maxiter'};
			for i = 1:numel(int_entries)
				if (isfield(settings,int_entries{i}))
					val = int32(getfield(settings,  int_entries{i})); %#ok<*GFLD>
					settings = setfield(settings, int_entries{i},val); %#ok<*SFLD>
				end
			end
			
			bool_entries = {'verbose','use_a_star','store_visit_time','store_parents' };
			for i = 1:numel(bool_entries)
				if (isfield(settings,bool_entries{i}))
					val = logical(getfield(settings,  bool_entries{i}));
					settings = setfield(settings, bool_entries{i},val);
				end
			end
		end

		function data = preprocess_data(self,data)
			data(self.disallowed_set) = inf;
		end

		% Create a mesh_map which is used internally to keep track on
		% start set, end_set and disallowed pixel/voxels.
		function self = create_mesh_map(self, start_set, end_set, disallowed)
			% Create a structure holding the start/end and allowed pixels
			% for the algorithm to visit.
			mesh_map = ones(self.problem_size,'uint8'); %#ok<*PROP>

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
				
				check = false(1,3);
				check(1) = info.pair < self.global_limit(1);
				check(2) = info.triplet < self.global_limit(2);
				check(3) = info.quad < self.global_limit(3);
				
				for i = 1:3
					if (self.global_limit(i) == inf)
						check(i) = true;
					end
				end

				feasible = all(check);
		end
	end

	methods

		% Constructor
		function self = Curve_extraction_base()
			addpath([fileparts(mfilename('fullpath')) filesep 'library']);
		end

		% Shortest path from the start set to the end set.
		function  [curve, cost, time, evaluations, visit_map ] = shortest_path(self)
			
			if (~any(self.mesh_map(:) == 3))
				error('The problem has no end set.');
			end

			settings = gather_settings(self);	
			data = self.preprocess_data(self.data);
			compile('curve_segmentation');

			[curve, ~, time, evaluations, visit_map] = ...
				curve_segmentation_mex(self.data_type, self.mesh_map, data, self.connectivity, settings);

			% Saving solution
			self.curve = curve;
			self.visit_map = visit_map;

			% When self.curve is set the cost is automatically updated with detailed info
			cost = self.cost;

			%% Supergradient
			% Check if curve is within limits otherwise run super-gradient optimization
			if (self.curve_feasibility_check(self.get_info() ))
				return;
			end

			% Variables
			% x(1): function value
			% x(2): two point penalty variable
			% x(3): three point penalty variable
			% x(4): four point penalty variable
			c = [-1; 0; 0; 0];
			A = zeros(0, 4);
			b = zeros(0, 1);
			lb = [-inf; self.range(1,1); self.range(2,1); self.range(3,1)];
			ub = [inf;  self.range(1,2); self.range(2,2); self.range(3,2)];

			projected_cost = nan;

			% Have we found any feasible solution?
			feasible_solution = false;

			% Determine which constraints are needed
			active_constraints = self.global_limit(1:3) < inf;

			initial_multipliers = settings.penalty(1:3);
			multipliers = initial_multipliers;

			limits = self.global_limit(1:3);
			limits(~active_constraints) = 0;

			for iter = 1:self.supergradient_iterations

				if iter > 1
					% Determine which point on the dual function to evaluate.
					options = optimset('Display','none');
					[x,~,exitflag] = linprog(c, A, b, [], [], lb, ub,[], options);
					
					if (exitflag < 0)
						error('err:constrained_linprog','Unable to solve a linear program inside constrained shortest path');
					end

					projected_cost = x(1);

					multipliers = x(2:4);
					multipliers(~active_constraints) = initial_multipliers(~active_constraints);

					settings.penalty(1:3) = multipliers;
		
					[sg.curve, ~, sg.time, sg.evaluations, sg.visit_map] = ...
						curve_segmentation_mex(self.data_type, self.mesh_map, self.data, self.connectivity, settings);

					time = time+sg.time;
					evaluations = evaluations + sg.evaluations;

					[cost,info] = self.curve_info(sg.curve, settings);
				else

					% Shortest path has already been run with initial penalty settings.
					cost = self.get_cost;
					info = self.get_info;
				end

				values(1,1) = info.pair;
				values(2,1) = info.triplet;
				values(3,1) = info.quad;

				% Constraint g(x) <= 0.
				g = values - limits;
				g(~active_constraints) = 0;

				% Dual function.
				% The solver returns  \min f(x) + \sum \multiplier * value
				% the dual function is defined as d(x) =  \min f(x) + \sum \multiplier*(value -limit)
				dual_val = cost.total - sum(multipliers.*limits);

				% Duality gap
				gap = projected_cost - dual_val;

				A = [A; 1, -g(:)']; %#ok<AGROW>
				b = [b; dual_val - sum(multipliers.*g)]; %#ok<AGROW>

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
					fprintf('Two point function: %g max: %g	penalty: %g \n', values(1), limits(1), multipliers(1));
					fprintf('Three point function: %g max: %g	penalty: %g. \n',  values(2), limits(2), multipliers(2));
					fprintf('Four point function: %g max: %g	penalty: %g. \n',  values(3), limits(3), multipliers(3));
					fprintf('---\n');
					drawnow();
				end

				if gap < self.supergradient_max_duality_gap
					break;
				end
			end

			if (~feasible_solution)
				warning('Unable to find a feasible solution. Returning optimal solution for the original cost function.');
			else
				if (self.override_penalty)
					self.penalty(1) = settings.penalty(1);
					self.penalty(2) = settings.penalty(2);
					self.penalty(3) = settings.penalty(3);
				end

				self.curve = curve;
				self.visit_map = visit_map;
			end
		end


		% Compute distance to every point from the start set
		function [distances, curve] = compute_all_distances(self)
      compile('curve_segmentation');

      data = self.preprocess_data(self.data);

			settings = gather_settings(self);
      settings.store_distances = true;
      settings.compute_all_distances = true;

			[curve, ~, ~, ~, ~, ~, distances] = ...
			 		 curve_segmentation_mex(self.data_type, self.mesh_map, data, self.connectivity, settings);

			if min(distances(self.disallowed_set))	 < 1e38
				error('Disallowed set should have infinite cost.');
			end

			distances(self.disallowed_set) = inf;

		end

		% Compute visitation_tree
		% Returns a matrix (named tree) where each entry contains the linear index to
		% the parent voxel when every pixel/voxel is visited.
    function  [tree, curve] = compute_visit_tree(self)

    	data = self.preprocess_data(self.data);

			settings = gather_settings(self);
      settings.store_parents = true;
      settings.compute_all_distances = true;
			
      compile('curve_segmentation');
			[curve, ~, ~, ~, ~, tree] = ...
			 		 curve_segmentation_mex(self.data_type, self.mesh_map, data, self.connectivity, settings);

			% If no curve is store this
			if isempty(self.curve)
				self.curve = curve;
			end
		end

		% Perform local optimization on a discrete shortest_path solution.
		% The user can perform local optimization on any curve
		% by setting obj.curve and then calling this method.
		function [curve,cost,time] = local_optimization(self)

			if isempty(self.curve)
				error('No curve stored, please set obj.curve or run obj.shortest_path()');
			end


			if (~any(self.end_set))
				error('The problem has no end set.');
			end

		curve = self.curve;

			% Add more points first time
			if ~self.checked_max_curve_segment_length
				curve = self.interpolate_more_points(curve, self.local_optimzation_max_curve_segment_length);
			end

			data = self.preprocess_data(self.data);

			settings = gather_settings(self);
			compile('local_optimization');

			[curve, ~, time, success] = ...
			local_optimization_mex(	self.data_type, data, curve, self.connectivity, settings);


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

			self.checked_max_curve_segment_length = true;
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
					new_curve = [new_curve; s+v*k]; %#ok<AGROW>
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
			assert(length(voxel_dimensions) > 1 && length(voxel_dimensions) < 5);
			assert(all(voxel_dimensions) > 0);

			if (length(voxel_dimensions) < 4)
				voxel_dimensions =[voxel_dimensions 1 1];
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

		function set.penalty(self, penalty)
			if ~(all(penalty) >= 0)
				error('Regularization coefficients must be non-negative')
			end

			self.penalty = penalty;
			self.reset_solution();
		end

		function set.power(self, power)
			if ~(all(power) >= 0)
				error('Regularization power must be non-negative')
			end

			self.power = power;
			self.reset_solution();
		end

		function set.local_limit(self, local_limit)
			if ~(all(local_limit) >= 0)
				error('Regularization local limit must be non-negative')
			end

			self.local_limit = local_limit;
			self.reset_solution();
		end

		function set.global_limit(self, global_limit)
			if ~(all(global_limit) > 0)
				error('Regularization local limit must be non-negative')
			end

			self.global_limit = global_limit;
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

		% Calculate curve cost on demand
		function cost = get_cost(self)
			if isempty(self.cached_cost)
				[self.cached_cost,self.cached_info] = self.curve_info(self.curve);
			end

			cost = self.cached_cost;
		end

		function info = get_info(self)
			if isempty(self.cached_cost)
				[self.cached_cost,self.cached_info] = self.curve_info(self.curve);
			end

			info = self.cached_info;
		end

		% Keeping the curve
		function reset_solution(self)
			self.cached_cost = [];
		end

		% Interface the start,end and disallowed set
		function allowed = get.disallowed_set(self)
			allowed = self.mesh_map == 0;
		end

		function set.disallowed_set(self, disallowed)
			assert(all(self.problem_size == size(disallowed))); %#ok<*MCSUP>

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
			self.checked_max_curve_segment_length = false;
			self.reset_solution();
		end

		function check_limit(~, limit)
			if (limit(1) < 0)
				error('Limit must be >= 0.');
			end
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
		% of the calculated with current three_points_power and four_points_power.
		function [cost,info] = curve_info(self, curve, settings)

			if nargin < 2
				curve = self.curve;
			end

			if (size(curve,2) == 2)
			   curve = [curve ones(size(curve,1),1)];
			end

			% Default
			cost.total = nan;
			cost.data = nan;
			cost.pair = nan;
			cost.triplet = nan;
			cost.quad = nan;
			
			info.data = nan;
			info.pair = nan;
			info.triplet = nan;
			info.quad = nan;

			
			if ~isempty(curve)
				if nargin < 3
					settings = gather_settings(self);
				end
				
				compile('curve_info');

				[cost.total, cost.data, cost.pair, cost.triplet, cost.quad, ...
				info.pair, info.triplet, info.quad] ...
				= curve_info_mex(self.data_type, self.data, curve,  ...
													self.connectivity, settings);
												
				info.data = cost.data;
			end
		end
	end
end