% Shortest path with length, curvature and torsions priors.
% ---
%
% Let
% f = data_cost(curve segment) 
% + length_penalty (length of curve segment)
% + curvature_penalty |curvature of curve segment|^curvature_power
% + torsion_penalty |torsion of curve segment|^torsion_power
%
% The data cost is given by nearest neighbor interpolation of the data matrix.
%
% This class minimizes:
% ----
% \int f ds
% s.t.
% \int (length of curve segment)ds <= length_global_limit
% \int (curvature_penalty | curvature of curve segment|^curvature_power)ds <= curvature_global_limit
% \int (torsion_penalty |torsion of curve segment|^torsion_power)ds <= torsion_global_limit.
%
classdef Curve_extraction < Curve_extraction_base
	
	properties (SetAccess = protected)
		cost;
		info;
		data_type = 'linear_interpolation';
	end
	
	properties
		% Virtual
		length_penalty;
		curvature_penalty;
		torsion_penalty;
		
		curvature_power;
		torsion_power;
		
		length_global_limit;
		curvature_global_limit;
		torsion_global_limit;
	end

	methods
		% Set, get methods
		
		% Penalty
		function set.length_penalty(self, length_penalty)
			self.penalty(1) = length_penalty;
		end
		
		function set.curvature_penalty(self, curvature_penalty)
			self.penalty(2) = curvature_penalty;
		end
		
		function set.torsion_penalty(self, torsion_penalty)
			self.penalty(3) = torsion_penalty;
		end
		
		function length_penalty = get.length_penalty(self)
			length_penalty = self.penalty(1);
		end
		
		function curvature_penalty = get.curvature_penalty(self)
			curvature_penalty = self.penalty(2);
		end
		
		function torsion_penalty = get.torsion_penalty(self)
			torsion_penalty = self.penalty(3);
		end
		
		% Power
		function set.curvature_power(self, curvature_power)
			self.power(2) = curvature_power;
		end
		
		function set.torsion_power(self, torsion_power)
			self.power(3) = torsion_power;
		end
		
		function curvature_power = get.curvature_power(self)
			curvature_power = self.power(2);
		end
		
		function torsion_power = get.torsion_power(self)
			torsion_power = self.power(3);
		end
				
		% Global limit
		function set.length_global_limit(self, length_global_limit)
			self.global_limit(1) = length_global_limit;
		end
		
		function set.curvature_global_limit(self, curvature_global_limit)
			self.global_limit(2) = curvature_global_limit;
		end
		
		function set.torsion_global_limit(self, torsion_global_limit)
			self.global_limit(3) = torsion_global_limit;
		end
		
		function length_global_limit = get.length_global_limit(self)
			length_global_limit = self.global_limit(1);
		end
		
		function curvature_global_limit = get.curvature_global_limit(self)
			curvature_global_limit = self.global_limit(2);
		end
		
		function torsion_global_limit = get.torsion_global_limit(self)
			torsion_global_limit = self.global_limit(3);
		end

		function cost = get.cost(self)
			base_cost = self.get_cost();
		
			cost.data = base_cost.data;
			cost.total = base_cost.total;
			cost.length = base_cost.pair;
			cost.curvature = base_cost.triplet;
			cost.torsion = base_cost.quad;
		end

		function info = get.info(self)
			base_info = self.get_info;
			
			info.data = base_info.data;
			info.length =  base_info.pair;
			info.curvature= base_info.triplet;
			info.torsion = base_info.quad;
		end
	end
	
	methods
		
		% Used when the data cost is defined by by performing linear interpolation a 2D or 3D matrix.
		function self = Curve_extraction(data, varargin)
			self = self@Curve_extraction_base;

			self.data = data;
			self.problem_size  = size(self.data);
			
			self.create_mesh_map(varargin{:});
			self.set_connectivity_by_radius(self.default_connectivity_radius);

			self.curvature_power = 2;
			self.torsion_power = 2;
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
			
			hold on;
			msgs = {};
			
			if (~isempty(curve))
				msgs{end+1} = sprintf('Cost; total: %g data: %g, length: %g curvature: %g, torsion %g.', ...
					self.cost.total, self.cost.data, self.cost.length, self.cost.curvature, self.cost.torsion);
			end
			
			msgs{end+1} = sprintf('Cost function; data + %g|length| + %g|curvature|^{%g} + %g|torsion|^{%g}.', ...
				self.length_penalty, self.curvature_penalty,self.curvature_power,self.torsion_penalty,self.torsion_power);
			
			if (~isempty(curve))
				msgs{end+1} = sprintf('Curve info: Length: %g |curvature|^{%g}: %g  |torsion|^{%g}: %g', ...
					self.info.length,  self.curvature_power, self.info.curvature,  self.torsion_power, self.info.torsion);
			end

			if (length(self.problem_size) == 2 && strcmp(self.data_type,'linear_interpolation'))
				cost_im = self.data;
				cost_im(self.mesh_map ~= 1) = -1;
				imagesc(double(cost_im))
				colormap gray(256); axis equal; axis off;

				axis ij;
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
	end
end