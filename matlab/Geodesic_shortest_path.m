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

			% if (~isempty(curve))
			% 	msgs{end+1} = sprintf('Cost; total: %g data: %g, length: %g curvature: %g, torsion %g.', ...
			% 		self.cost.total, self.cost.data, self.cost.length, self.cost.curvature, self.cost.torsion);
			% end

			% msgs{end+1} = sprintf('Cost function; data + %g|length| + %g|curvature|^{%g} + %g|torsion|^{%g}.', ...
			% 		self.length_penalty, self.curvature_penalty, self.curvature_power, self.torsion_penalty, self.torsion_power);

			% if (~isempty(curve))
			% 	msgs{end+1} = sprintf('Curve info: Length: %g |curvature|^{%g}: %g  |torsion|^{%g}: %g', ...
			% 		self.info.length,  self.curvature_power, self.info.curvature,  self.torsion_power, self.info.torsion);
			% end

			cost_im = self.data;
			cost_im(self.mesh_map ~= 1) = -1;
			imagesc(double(cost_im))
			colormap gray(256); axis equal; axis off;
			axis ij;


			if (~isempty(curve))
				title({msgs{:}});
				cmap = jet(3);

				plot(curve(:,2),curve(:,1),'r-' , 'linewidth',2)
				plot(curve(1,2),curve(1,1),'ko','MarkerFaceColor', cmap(2,:), 'MarkerSize',5);
				plot(curve(end,2),curve(end,1),'ko','MarkerFaceColor', cmap(3,:), 'MarkerSize',5);

				legend('Curve','Start','End','Location', 'EastOutside');
			else
				fprintf('No solution stored, please run obj.shortest_path() \n');
			end
		end
	end
end