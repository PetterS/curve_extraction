% Shortest path with length, curvature and torsions priors.
% ---
%
% Let
% f = data_cost(curve segment) 
% + length_penalty (length of curve segment)
% + curvature_penalty |curvature of curve segment|^curvature_power
% + torsion_penalty |torsion of curve segment|^torsion_power
%
% The data cost is explicitly defined for each direction by the data matrix (see an example in /examples)
%
% This class minimizes:
% ----
% \int f ds
% s.t.
% \int (length of curve segment)ds <= length_global_limit
% \int (curvature_penalty | curvature of curve segment|^curvature_power)ds <= curvature_global_limit
% \int (torsion_penalty |torsion of curve segment|^torsion_power)ds <= torsion_global_limit.
%
classdef Curve_extraction_edge < Curve_extraction
	methods
		
		function self = Curve_extraction_edge(data, connectivity, varargin)
			self = self@Curve_extraction(1,true);

			sz = size(data);
			self.problem_size = sz(1:end-1);
			
			self.data = data;
			self.connectivity = int32(connectivity);
			self.data_type = 'edge';

			self.create_mesh_map(varargin{:});
		end
		
		function local_optimization(~)
			error('Local optimization is not possible when the data cost is given explicitly');
		end

		function plot_curve(self,curve)

			if nargin < 2
				curve = self.curve;
			end

			% Just display start, end, allowed set.
			if (length(self.problem_size) == 2)
				imagesc(3-self.mesh_map);
				colormap(gray(4)); axis equal; axis off;
				axis ij;
			end

			plot_curve@Curve_extraction(self,curve);
		end
	end

	methods (Access = protected)
		function data = preprocess_data(self,data)
	
			if ~any(self.disallowed_set(:))
				return
			end
			
			if (numel(self.problem_size) == 2)
				for d = 1:size(self.data,3)
					slice = data(:,:,d);
					slice(self.disallowed_set) = inf;
					data(:,:,d) = slice;
				end
			end
	
			if (numel(self.problem_size) == 3)
				for d = 1:size(self.data,3)
					slice = data(:,:,:,d);
					slice(self.disallowed_set) = inf;
					data(:,:,:,d) = slice;
				end
			end
		end
	end
end