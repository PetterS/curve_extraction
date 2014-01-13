% Class handling curve extraction
% Johannes
classdef Curve_extraction < handle
	
	% Settings
	properties
		length_penalty = 0;
		curvature_penalty = 0,
		torsion_penalty = 0,
		power_curvature = 2.0;
		power_torsion = 2.0;
		use_a_star = true;
		verbose = false;
		store_visit_time = false;
		num_threads = int32(1);

		unary =  [];

		% Created in the constructor.
		connectivity = []; 

		% Local optimization
		maxiter = 1000;
		descent_method = 'lbfgs';
		
		voxel_dimensions = [1 1 1];
		
		curve = [];
	end
	
	% Stored by the solver
	properties (SetAccess = protected)
		unary_type = '';
		time = nan;
		cost = nan;
		evaluations = nan,
		mesh_map = [];
		
		% Note: For length regularization the visit map runs from end set to start set
		% this is because the same code get lower bound for A*.
		visit_map = [];
	end
	
	properties (Hidden)
		problem_size;
	end

	methods (Access = protected)
		% Convert settings to format supported by the mex file.
		function settings = gather_settings(self)
			settings.length_penalty = self.length_penalty;
			settings.curvature_penalty = self.curvature_penalty;
			settings.torsion_penalty = self.torsion_penalty;
			settings.power_curvature = self.power_curvature;
			settings.use_a_star = self.use_a_star;
			settings.verbose = self.verbose;
			settings.unary_type = self.unary_type;
			settings.maxiter = self.maxiter;
			settings.store_visit_time =  self.store_visit_time;
			settings.descent_method = self.descent_method;	
			settings.voxel_dimensions = self.voxel_dimensions;
			settings.num_threads = self.num_threads;
		end

		function interpolation_constructor(self, unary, varargin)		
			self.problem_size = size(unary);
			self.unary = unary;

			default_radius = 4;
			self.set_connectivity_by_radius(default_radius);

			self.create_mesh_map(varargin{:});
		end

		function edge_constructor(self, unary, connectivity, varargin)		
			self.connectivity = int32(connectivity);
			self.problem_size = size(unary);
			self.unary = unary;

			sz = size(unary);
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
		function self = Curve_extraction(unary_type, varargin)
			addpath([fileparts(mfilename('fullpath')) filesep 'library']);
			
			self.unary_type = unary_type;

			% Send the input to the correct constructor:
			switch unary_type
				case 'linear_interpolation'
					self.interpolation_constructor(varargin{:});
				case 'edge'
					self.edge_constructor(varargin{:});
			end

			if (length(self.problem_size) ~= 2) && (length(self.problem_size) ~= 3)
				error('Only two- and three-dimensional problems supported.');
			end	

		end
		
		% Solution cost decomposed in the different terms.
		function cost = curve_info(self)

			if ~strcmp(self.unary_type,'linear_interpolation')
				error('curve_info is only supported linear_interpolation unary costs _At the moment_');
			end

			settings = gather_settings(self);
			cost = curve_info(self.unary, self.curve, settings);
		end
		
				
		function  [curve, cost, time, evaluations, visit_map ] = solve(self)
			
			if ~strcmp(self.unary_type,'linear_interpolation')
				error('Solve is only supported for linear_interpolation unary costs _At the moment_');
			end

			settings = gather_settings(self);
						
			[curve, cost, time, evaluations, visit_map] = ...
			 		 curve_segmentation(self.mesh_map, self.unary, self.connectivity, settings);
			
			% Saving solution
			self.curve = curve;
			self.time = time;
			self.cost =  cost;
			self.evaluations = evaluations;
			self.visit_map = visit_map;
        end
        
    function  [tree] = compute_tree(self)
			settings = gather_settings(self);
      settings.store_parents = true;
			
			[~, ~, ~, ~, ~, tree] = ...
			 		 curve_segmentation(self.mesh_map, self.unary, self.connectivity, settings);
		end
		
		% Move away from discretized solution and find a local optimum.
		function [curve,cost,time] = local_optimization(self)
			
			if ~strcmp(self.unary_type,'linear_interpolation')
				error('Local optimization is only supported  for linear_interpolation unary costs');
			end

			if isempty(self.curve)
				fprintf('No curve stored, running the solver \n');
				self.solve()
			end
			
			
			settings = gather_settings(self);
			[curve, cost, time] = local_optimization(self.mesh_map, self.unary, self.curve, settings);
			
			if (cost > self.cost)
				warning('Unable to find a better local optima');
				curve = self.curve;
				cost = self.cost;
			else
				% Saving solution
				self.curve = curve;
				self.cost = cost;
			end
		end
		
		% Draw current solution (only supports 2D curves)
		function display(self)
			details(self);

			if (~ismatrix(self.unary))
				if (~isempty(self.curve))
					plot3(self.curve(:,1), self.curve(:,2), self.curve(:,3),'-r');
					title(sprintf('Solution cost: %g \n', self.cost));
					fprintf('Solution cost: %g \n', self.cost);
				end
			else
				
				figure();
				cost_im = self.unary;
				cost_im(self.mesh_map ~= 1) = -1;
				imagesc(double(cost_im))
				colormap gray(256); hold on; axis equal; axis off;
				
				if (~isempty(self.curve))
					plot(self.curve(:,2),self.curve(:,1),'-r' , 'linewidth',5)
					title(sprintf('Solution cost: %g \n', self.cost));
					fprintf('Solution cost: %g \n', self.cost);
				else
					fprintf('No solution stored, please run obj.solve() \n');
				end								
			end
		end
		
		%% Set functions
		function set_connectivity_by_radius(self, radius)

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
		
		function set.connectivity(self, connectivity)

			% Padding
			if (size(connectivity,2) == 2)
				connectivity = [connectivity zeros(size(connectivity,1),1)];
			end
			
			assert(size(connectivity,2) == 3);

			self.connectivity = connectivity;
		end

		function set.unary(self, unary)
			self.unary = unary;
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
		
		% Keeping the cuvre
		function reset_solution(self)
			self.cost = nan;
		end
		
		function set.unary_type(self, unary_type)

			switch unary_type
				case 'linear_interpolation'
					%ok
				case 'edge'
					%ok
				otherwise
					error(sprintf('Possible unary types = {linear_interpolation, edge}'));
			end
			
			self.unary_type = unary_type;
		end
		
		% Interface the start,end and disallowed set
		function allowed = get_disallowed_set(self)
			allowed = self.mesh_map == 0;
		end
		
		function set_disallowed_set(self, disallowed)
			assert(all(self.problem_size == size(disallowed)));
			
			self.mesh_map(~disallowed & self.mesh_map==0) = 1;
			self.mesh_map(disallowed) = 0;
		end
		
		function start_set = get_start_set(self)
			start_set = self.mesh_map == 2;
		end
	
		function set_start_set(self, start_set)
			assert(all(self.problem_size == size(start_set)));
			self.mesh_map(start_set) = 2;
		end
		
		function end_set = get_end_set(self)
			end_set = self.mesh_map == 3;
		end
		
		function set_end_set(self, end_set)
			assert(all(self.problem_size == size(end_set)));
			self.mesh_map(end_set) = 3;
		end
		
		function set.mesh_map(self, mesh_map)
			
			if (~any(mesh_map(:) == 2))
				error('The problem has no start set.');
			end

			if (~any(mesh_map(:) == 3))
				error('The problem has no end set.');
			end
			
			self.mesh_map = mesh_map;
		end
		
		function set.curve(self, curve)
			
			assert(size(curve,2) == length(self.problem_size));
			self.curve = curve;
		end
			
	end
end