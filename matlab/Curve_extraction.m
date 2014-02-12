% Class handling curve extraction
% Johannes
classdef Curve_extraction < handle

	% Settings
	properties
		length_penalty = 0;
		curvature_penalty = 0,
		torsion_penalty = 0,
		curvature_power = 2.0;
		torsion_power = 2.0;
		use_a_star = true;
		verbose = false;
		store_visit_time = false;
		num_threads = int32(min(2,feature('numThreads')));

		data =  [];

		% Created in the constructor.
		connectivity = [];

		% Local optimization
		maxiter = 1000;
		descent_method = 'lbfgs';

		voxel_dimensions = [1 1 1];
		curve = [];

		% "Virtual"
		start_set;
		end_set;
		disallowed_set;
	end

	% Stored by the solver
	properties (SetAccess = protected)
		data_type = '';
		time = nan;
		cost;
		info;
		evaluations = nan,

		% Note: For length regularization the visit map runs from end set to start set
		% this is because the same code get lower bound for A*.
		visit_map = [];
	end

	properties (Hidden)
		problem_size;
		cached_cost;
		cached_info;
		mesh_map = [];
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
			settings.maxiter = self.maxiter;
			settings.store_visit_time =  self.store_visit_time;
			settings.descent_method = self.descent_method;
			settings.voxel_dimensions = self.voxel_dimensions;
			settings.num_threads = self.num_threads;
		end

		function geodesic_constructor(self, data, varargin)
	
			% Default length penalty is ~= 0 for this data type
			self.length_penalty = 1;

			if ~ismatrix(data)
				error('Only 2D dimensional data supported by geodesic shortest path');
			end
			
			self.interpolation_constructor(data, varargin{:});
		end
		
		function interpolation_constructor(self, data, varargin)
			self.problem_size = size(data);
			self.data = data;

			default_radius = 3;
			self.set_connectivity_by_radius(default_radius);

			self.create_mesh_map(varargin{:});
		end

		function edge_constructor(self, data, connectivity, varargin)
			self.connectivity = int32(connectivity);
			self.problem_size = size(data);
			self.data = data;

			sz = size(data);
			self.problem_size = sz(1:end-1);

			self.create_mesh_map(varargin{:});
		end


		function self = create_mesh_map(self, start_set, end_set, disallowed)
			% Create a structure holding the start/end and allowed pixels
			% for the algorithm to visit.
			mesh_map = ones(self.problem_size,'uint8');

			if nargin > 3
				assert(all(self.problem_size == size(disallowed)));
				mesh_map(disallowed) = 0;
			end

			% Check input
			assert(all(self.problem_size == size(start_set)));
			assert(all(self.problem_size == size(end_set)));

			if (any(start_set(:) & end_set(:)))
				error('Some voxels are both in the start and end set');
			end

			mesh_map(start_set) = 2;
			mesh_map(end_set) = 3;

			% Save
			self.mesh_map = mesh_map;
		end
	end
	methods
		% Find global solution in the mesh implicitly defined by the connectivity.
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

		function  [curve, cost, time, evaluations, visit_map ] = solve(self)
			settings = gather_settings(self);

			[curve, total_cost, time, evaluations, visit_map] = ...
			 		 curve_segmentation(self.data_type, self.mesh_map, self.data, self.connectivity, settings);

			% Saving solution
			self.curve = curve;
			self.time = time;
			self.evaluations = evaluations;
			self.visit_map = visit_map;

			% When self.curve is set the cost is automatically updated
			% with detailed info
			cost = self.cost;
		end

    function  [tree, curve] = compute_tree(self)
			settings = gather_settings(self);
      settings.store_parents = true;

			[curve, ~, ~, ~, ~, tree] = ...
			 		 curve_segmentation(self.mesh_map, self.data, self.connectivity, settings);

			% If no curve is store this
			if isempty(self.curve)
				self.curve = curve;
			end
		end

		% Move away from discretized solution and find a local optimum.
		function [curve,cost,time] = local_optimization(self)

			if ~strcmp(self.data_type,'linear_interpolation')
				error('Local optimization is only supported  for linear_interpolation data costs');
			end

			if isempty(self.curve)
				fprintf('No curve stored, running the solver \n');
				self.solve()
			end

			settings = gather_settings(self);
			[curve, ~, ~, success] = local_optimization(self.mesh_map, self.data, self.curve, settings);

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

		% Draw current solution (only supports 2D curves)
		function display(self)
			details(self);
			clf; hold on;

			msgs = {};
			
			if (~isempty(self.curve))
				msgs{end+1} = sprintf('Cost; total: %g data: %g, length: %g curvature: %g, torsion %g.', ...
					self.cost.total, self.cost.data, self.cost.length, self.cost.curvature, self.cost.torsion);
			end
			
			msgs{end+1} = sprintf('Penalty; %g|length| + %g|curvature|^{%g} + %g|torsion|^{%g}.', ...
					self.length_penalty, self.curvature_penalty, self.curvature_power, self.torsion_penalty, self.torsion_power);

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

			if (~isempty(self.curve))
				title({msgs{:}});

				cmap = jet(3);

				% Draw the stored solution.
				if (length(self.problem_size) == 3)
					plot3(self.curve(:,1), self.curve(:,2), self.curve(:,3),'-r', 'linewidth', 2);

					plot3(self.curve(1,1),self.curve(1,2),self.curve(1,3),'ko', 'MarkerFaceColor',cmap(2,:), 'MarkerSize',5);
					plot3(self.curve(end,1),self.curve(end,2),self.curve(end,3),'ko','MarkerFaceColor', cmap(3,:), 'MarkerSize',5);

					legend('Curve','Start','End','Location', 'EastOutside');
					view(3)
				else

					plot(self.curve(:,2),self.curve(:,1),'r-' , 'linewidth',2)
					plot(self.curve(1,2),self.curve(1,1),'ko','MarkerFaceColor', cmap(2,:), 'MarkerSize',5);
					plot(self.curve(end,2),self.curve(end,1),'ko','MarkerFaceColor', cmap(3,:), 'MarkerSize',5);

					legend('Curve','Start','End','Location', 'EastOutside');
				end
			else
				fprintf('No solution stored, please run obj.solve() \n');
			end
		end

		%% Set functions
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
			assert(length(voxel_dimensions) == 2 || length(voxel_dimensions == 3));
			assert(all(voxel_dimensions) > 0);

			if (length(voxel_dimensions) == 2)
				voxel_dimensions =[voxel_dimensions 1];
			end

			self.voxel_dimensions = voxel_dimensions;
			self.reset_solution();
		end


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

		% Connectivity
		% Each row defines offset from a voxel,
		% e.g: For 2D data connectivity = [1 0; 0 1; -1 0; 0 -1],
		% defines standard 4-connectivity.
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

		function set.maxiter(self, maxiter)
			if (maxiter < 0)
				error('Number of iterations for local optimization must be positive.');
			end

			self.maxiter = maxiter;
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

			if (~any(mesh_map(:) == 3))
				error('The problem has no end set.');
			end

			self.mesh_map = mesh_map;
			self.reset_solution();
		end

		function set.curve(self, curve)
			assert(size(curve,2) == length(self.problem_size));
			self.curve = curve;
			self.reset_solution();
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