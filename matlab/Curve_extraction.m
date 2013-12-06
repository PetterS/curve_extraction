% Class handling curve extraction
% Johannes
classdef Curve_extraction < handle
	properties
		info = [];
		curve = [];
		cost = 0;
		length_penalty = 0;
		curvature_penalty = 0,
		torsion_penalty = 0,
		power_curvature = 2.0;
		power_torsion = 2.0;
		regularization_radius = 4.0
		use_a_star = true;
		VERBOSE = false;
		store_visit_time = false;
		unary_type = 'linear';
		unary =  [];
		maxiter = 100;
		mesh_map = [];
	end
	
	methods (Access = protected)
		% Convert settings to format supported by the mex file.
		function settings = gather_settings(self)
			settings.length_penalty = self.length_penalty;
			settings.curvature_penalty = self.curvature_penalty;
			settings.torsion_penalty = self.torsion_penalty;
			settings.power_curvature = self.power_curvature;
			settings.regularization_radius = self.regularization_radius;
			settings.use_a_star = self.use_a_star;
			settings.VERBOSE = self.VERBOSE;
			settings.unary_type = self.unary_type;
			settings.maxiter = self.maxiter;
			settings.store_visit_time =  self.store_visit_time;
		end
	end
	methods
		% Find global solution in the mesh impicitly defined by the connectivity.
		function self = Curve_extraction(mesh_map, unary)
			addpath([fileparts(mfilename('fullpath')) filesep 'library']);
            
			self.mesh_map = mesh_map;
			self.unary = unary;
		end
		
		% Solution cost decomposed in the diffrent terms.
		function cost = curve_info(self)
			settings = gather_settings(self);
			cost = curve_info(self.unary, self.curve, settings);
		end
		
		function solution = solve(self)
			settings = gather_settings(self);

			solution = curve_segmentation(self.mesh_map, self.unary, settings);
			
			% Saving solution
			self.curve = solution.path;
			self.info = rmfield(solution,'path');
			self.cost = self.info.cost;
		end
		
		% Move away from discretized solution and find a local optimum.
		function curve = local_optimization(self)
			
			if isempty(self.curve)
				fprintf('No curve stored, running the solver \n');
				self.solve()
			end
			
			settings = gather_settings(self);
			[curve, info] = local_optimization(self.mesh_map, self.unary, self.curve, settings);

			% Saving solution
			self.curve = curve;
			self.cost = info.cost;
		end
		
		% Draw current solution (only supports 2D curves)
		function display(self)
			if (~ismatrix(self.unary))
				fprintf('This function can only draw solution to 2D problems');
			else
				
				figure();
				cost_im = self.unary;
				cost_im(self.mesh_map ~= 1) = -1;
				imagesc(double(cost_im))
				colormap gray(256); hold on; axis equal; axis off;
				
				if (~isempty(self.curve))
					plot(self.curve(:,2),self.curve(:,1),'-r' , 'linewidth',5)
					title(sprintf('Solution cost: %g \n', self.cost));
				else
					fprintf('No solution stored, please run obj.solve() \n');
				end
			end
		end
		
		%% Set functions
		function set.mesh_map(self, mesh_map)
			self.mesh_map = int32(mesh_map);
			self.erase_solution();
		end
		
		function set.unary(self, unary)
			self.unary = unary;
			self.erase_solution();
		end
		
		function set.length_penalty(self, length_penalty)
			if (length_penalty < 0)
				error('Regularization coefficents must non-negative');
			end
			
			self.length_penalty = length_penalty;
			self.erase_solution();
		end
		
		function set.maxiter(self, maxiter)
			if (maxiter < 0)
				error('Number of iterations for local optimization must be positive.');
			end
			
			self.maxiter = maxiter;
		end
		
		function set.curvature_penalty(self, curvature_penalty)
			if (curvature_penalty < 0)
				error('Regularization coefficents must non-negative');
			end
			
			self.curvature_penalty = curvature_penalty;
			self.erase_solution();
		end
		
		function set.VERBOSE(self, VERBOSE)
			if ~isa(VERBOSE,'logical')
				error('VERBOSE must either be true or false');
			end
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
		end
		
		function set.torsion_penalty(self, torsion_penalty)
			if (torsion_penalty < 0)
				error('Regularization coefficents must non-negative');
			end
			
			self.torsion_penalty = torsion_penalty;
			self.erase_solution();
		end
		
		function erase_solution(self)
			self.info = [];
			self.curve = [];
		end
		
		function set.unary_type(self, unary_type)
			if ~strcmp(unary_type,'linear')
				error('Currently only linear dataterm supported.');
			end
			
			self.unary_type = unary_type;
		end
	end
end