% Geodesic shortest path
%
% This class find the geodesic shortest path between two sets on a surface.
% The surface is defined as the function graph of discrete depth map.
%
classdef Geodesic_shortest_path < Curve_extraction_base
	
	properties (SetAccess = protected)
		cost;
		info;
		data_type = 'geodesic';
	end

	methods
		function cost = get.cost(self)
			base_cost = self.get_cost();
			cost = base_cost.pair;
			
		end

		function info = get.info(self)
			base_info = self.get_info();
			info = base_info.pair;
		end
	end

	methods
		function self = Geodesic_shortest_path(data, varargin)
			self = self@Curve_extraction_base;

			if ~ismatrix(data)
				error('Only 2D dimensional data supported by geodesic shortest path');
			end

			self.data = data;
			self.problem_size  = size(self.data);
			self.create_mesh_map(varargin{:});
			self.set_connectivity_by_radius(self.default_connectivity_radius);

			self.penalty(1) = 1;
			
			% No caching, parallelization is effective.
			self.num_threads = int32(feature('numThreads'));
		end

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


			cost_im = self.data;
			cost_im(self.mesh_map ~= 1) = -1;
			imagesc(double(cost_im))
			colormap gray(256); axis equal; axis off;
			axis ij;


			if (~isempty(curve))
				title({msgs{:}});
				cmap = jet(3);

				title(sprintf('Curve length: %g.', self.length()));

				plot(curve(:,2),curve(:,1),'r-' , 'linewidth',2)
				plot(curve(1,2),curve(1,1),'ko','MarkerFaceColor', cmap(2,:), 'MarkerSize',5);
				plot(curve(end,2),curve(end,1),'ko','MarkerFaceColor', cmap(3,:), 'MarkerSize',5);

				legend('Curve','Start','End','Location', 'EastOutside');
			else
				fprintf('No solution stored, please run obj.shortest_path() \n');
			end
		end

		function plot_3d(self, style)
			% Add more points first time
			if ~self.checked_max_curve_segment_length
				self.curve = self.interpolate_more_points(self.curve, self.local_optimzation_max_curve_segment_length);
			end

			if isempty(self.curve)
				fprintf('No solution stored, please run obj.shortest_path() \n');
				return;
			end

			title(sprintf('Curve length: %g.', self.length()));
			curve3d = self.curve;

			for ind = 1:size(self.curve, 1)
				x = self.curve(ind, 1);
				y = self.curve(ind, 2);
				z = self.voxel_dimensions(3) * interp2(self.data, x, y, 'linear');
				curve3d(ind, 3) = z;
			end

			cmap = jet(3);
			plot3(curve3d(:,1), curve3d(:,2), curve3d(:,3), style , 'linewidth',2)
			plot3(curve3d(1,1), curve3d(1,2), curve3d(1,3), 'ko','MarkerFaceColor', cmap(2,:), 'MarkerSize',5);
			plot3(curve3d(end,1), curve3d(end,2), curve3d(end,3), 'ko','MarkerFaceColor', cmap(3,:), 'MarkerSize',5);

			legend('Curve','Start','End','Location', 'EastOutside');
		end

		function L = length(self)
			% Add more points first time
			if ~self.checked_max_curve_segment_length
				self.curve = self.interpolate_more_points(self.curve, self.local_optimzation_max_curve_segment_length);
			end

			L = self.info;
		end
	end
end