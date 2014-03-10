% Unit tests for the MATLAB code
%
% Perform all test by
% run(unit_tests)

classdef unit_tests < matlab.unittest.TestCase
	
	properties
		linear_obj;
		linear_obj_large;
		edge_obj;
		
		rng_seed = 0;
		problem_size = [9 8 5];
		problem_size_large = [50 50 25];
	end
	
	methods
		function switch_start_and_end_set(~, curve_obj)
			tmp = curve_obj.start_set;
			curve_obj.start_set = curve_obj.end_set;
			curve_obj.end_set = tmp;
		end

		% Computes the visit tree for C and checks that it is sane.
		function test_visit_tree(obj, C)
			[tree, ~] = C.compute_visit_tree();
			obj.assertEqual(size(tree), obj.problem_size);

			% Verify that every node in the tree points to something valid.
			for ind = 1:numel(tree)
				obj.assertGreaterThanOrEqual(tree(ind), 0);
				obj.assertLessThanOrEqual(tree(ind), prod(obj.problem_size));
			end

			% Verify that it is actually a tree.
			check = -ones(size(tree));
			for ind = 1:numel(tree)
				% Follow this path to the root.
				ind2 = ind;
				while tree(ind2) ~= 0 % The root node.
					% Next indices
					[i2, j2, k2] = ind2sub(size(tree), tree(ind2));
					% We cannot loop back.
					obj.assertLessThan(check(i2, j2, k2), ind);
					% Store that we are visiting this node.
					check(i2, j2, k2) = ind;
					% Go to the node.
					ind2 = sub2ind(size(tree), i2, j2, k2);
				end
			end
		end
	end
	
	% Setup
	methods (TestMethodSetup)
		function create_large_linear_obj(testCase)
			
			rng(testCase.rng_seed);
			data = rand(testCase.problem_size_large);
			
			start_set = false(testCase.problem_size_large);
			end_set = false(testCase.problem_size_large); 
			
			start_set(:,:,1) = true;
			end_set(:,:,end) = true;
			
			data_type = 'linear_interpolation';
			testCase.linear_obj_large = Curve_extraction(data_type,data, start_set, end_set);
			
		end
		
		function create_linear_obj(testCase)
			
			% Generate random data
			rng(testCase.rng_seed);
			data = rand(testCase.problem_size);
			
			start_set = false(testCase.problem_size);
			end_set = false(testCase.problem_size);
			disallowed_set = false(testCase.problem_size);
						
			start_set(5:6,5:6,1) = true;
			end_set(5:6,5:6,end) = true;
			disallowed_set(rand(testCase.problem_size) < 0.2) = 1;		

			data_type = 'linear_interpolation';
			testCase.linear_obj = Curve_extraction(data_type,data, start_set, end_set, disallowed_set);
		end
	end
	
	methods (TestMethodTeardown)
		function delete_objs(obj)
			clear obj.linear_obj;
			clear obj.linear_obj_large;
		end
	end
	
	methods (Test)
		
		% Test different regularizations
		function linear_data_cost(obj)
			tol = 1e-4;
			C = obj.linear_obj;
			
			%% Length
			C.set_connectivity_by_radius(3);
			C.length_penalty = 0.05;
			
			[curve,cost] = C.shortest_path();
			obj.verifyEqual(cost.total, 1.146388956723672, 'AbsTol', tol);
			obj.verifyEqual(cost.data, 0.889075738226573, 'AbsTol', tol);
			obj.verifyEqual(cost.length,  0.257313218497099, 'AbsTol', tol);
			obj.verifyEqual(cost.curvature, 0, 'AbsTol', tol);
			obj.verifyEqual(cost.torsion, 0, 'AbsTol', tol);
			obj.verifyEqual(size(curve,1), 5);


			%% |curvature|^2
			C.curvature_power = 2;
			C.curvature_penalty = 0.05;
			[curve,cost] = C.shortest_path();

			
			obj.verifyEqual(cost.total, 1.243813238105167, 'AbsTol', tol);
			obj.verifyEqual(cost.data,  0.909633257339023, 'AbsTol', tol);
			obj.verifyEqual(cost.length,  0.243185165257814, 'AbsTol', tol);
			obj.verifyEqual(cost.curvature,   0.090994815508330, 'AbsTol', tol);
			obj.verifyEqual(cost.torsion, 0, 'AbsTol', tol);
			obj.verifyEqual(size(curve,1), 4);
			
			%% |torsion|^2
			C.torsion_penalty = 1e-4;
			[curve,cost] = C.shortest_path();
			
			obj.verifyEqual(cost.total,  1.245438977306692, 'AbsTol', tol);
			obj.verifyEqual(cost.data,   0.909633257339023, 'AbsTol', tol);
			obj.verifyEqual(cost.length,   0.243185165257814, 'AbsTol', tol);
			obj.verifyEqual(cost.curvature,  0.090994815508330, 'AbsTol', tol);
			obj.verifyEqual(cost.torsion, 0.001625739201525, 'AbsTol', tol);
			obj.verifyEqual(size(curve,1), 4);			
		end
		
		%% A*
		% Using A* should always give the same solution as running the code
		% without it.
		function a_star(obj)
			C = obj.linear_obj;
			C.torsion_penalty = 0;
			tol = 1e-10;
			
			% Curvature
			for length_penalty = [0 0.1 1.38 4.23]
				for curvature_penalty = [0 0.2 0.33 0.46];
					C.length_penalty = length_penalty;
					C.curvature_penalty = curvature_penalty;
					
					C.use_a_star = false;
					[curve,cost] = C.shortest_path;
					
					C.use_a_star = true;
					[a_curve, a_cost] = C.shortest_path;
					
					obj.verifyEqual(numel(curve),numel(a_curve));
					obj.verifyEqual(cost.total, a_cost.total, 'AbsTol', tol);
				end
			end
			
			% Torsion
			C.length_penalty = 0.28;
			C.curvature_penalty = 0.43;
			for torsion_penalty = [0 0.2 0.33 0.46];

					C.torsion_penalty = torsion_penalty;
					C.use_a_star = false;
					[curve,cost] = C.shortest_path;
					
					C.use_a_star = true;
					[a_curve, a_cost] = C.shortest_path;
					
					obj.verifyEqual(numel(curve),numel(a_curve));
					obj.verifyEqual(cost.total, a_cost.total, 'AbsTol', tol);
			end
		end
		
		%% eps penalty 
		% Adding very small penalty to curvature or torsion will call a different solver
		% but should not change the resulting curve.
		function eps_penalty_increase(obj)
			C = obj.linear_obj;
			C.set_connectivity_by_radius(2.5);
			C.length_penalty = 1.1;
			C.curvature_penalty = 0;
			C.torsion_penalty = 0;
			tol = 1e-4;
			[curve,cost] = C.shortest_path();
			
			% Curvature
			C.curvature_penalty = 1e-100;
			[curve2,cost2] = C.shortest_path();
			obj.verifyEqual(numel(curve),numel(curve2));
			obj.verifyEqual(cost.total, cost2.total, 'AbsTol', tol);
			
			% Torsion
			C.curvature_penalty = 0;
			C.torsion_penalty = 1e-100;
			[curve3,cost3] = C.shortest_path();
			obj.verifyEqual(numel(curve),numel(curve3));
			obj.verifyEqual(cost.total, cost3.total, 'AbsTol', tol);
		end
		
		% Switching start and end set
		% should result in the same (but inverted)  curve
		function symmetric(obj)
			C = obj.linear_obj;
			tol = 1e-8;
			
			C.length_penalty = 0.48;
			C.curvature_penalty = 0;
			C.torsion_penalty = 0;
			
			%% Length
			[curve,cost] = C.shortest_path();
			obj.switch_start_and_end_set(C);
			[s_curve,s_cost] = C.shortest_path();
			obj.verifyTrue( all(all( curve == flipdim(s_curve,1) )) );
			obj.verifyEqual(cost.total, s_cost.total, 'AbsTol', tol);
			
			%% Curvature
			C.curvature_penalty = 0.28;
			[curve,cost] = C.shortest_path();
			obj.switch_start_and_end_set(C);
			[s_curve,s_cost] = C.shortest_path();
			obj.verifyTrue( all(all( curve == flipdim(s_curve,1) )) );
			obj.verifyEqual(cost.total, s_cost.total, 'AbsTol', tol);
			
			%% Torsion
			C.torsion_penalty = 0.11;
			[curve,cost] = C.shortest_path();
			obj.switch_start_and_end_set(C);
			[s_curve,s_cost] = C.shortest_path();
			obj.verifyTrue( all(all( curve == flipdim(s_curve,1) )) );
			obj.verifyEqual(cost.total, s_cost.total, 'AbsTol', tol);
			
			obj.switch_start_and_end_set(C);
		end
		
		%% Increase regularization cannot lead to lower cost
		function increased_regularization(obj)
			C = obj.linear_obj;
			C.set_connectivity_by_radius(3);
			C.length_penalty = 0;
			C.curvature_penalty = 0;
			C.torsion_penalty = 0;
			
			import matlab.unittest.constraints.IsGreaterThanOrEqualTo;
			
			%% Length
			prev_cost =0;
			for pen = linspace(0,10,10)
				C.length_penalty = pen;
				[~,cost] =  C.shortest_path;
				obj.verifyThat(cost.total, IsGreaterThanOrEqualTo(prev_cost));
				prev_cost = cost.total;
			end
			
			%% Curvature
			prev_cost =0;
			C.length_penalty = 0;
			for pen = linspace(0,10,10)
				C.curvature_penalty = pen;
				[~,cost] =  C.shortest_path;
				obj.verifyThat(cost.total, IsGreaterThanOrEqualTo(prev_cost));
				prev_cost = cost.total;
			end
			
			
			%% Torsion
			prev_cost =0;
			C.curvature_penalty = 0;
			for pen = linspace(0,10,10)
				C.torsion_penalty = pen;
				[~,cost] =  C.shortest_path;
				obj.verifyThat(cost.total, IsGreaterThanOrEqualTo(prev_cost));
				prev_cost = cost.total;
			end
		end
		
		%% Compute all distances tests
		function all_distances(obj)
			C = obj.linear_obj;
			distances = C.compute_all_distances;
			obj.verifyEqual(sum(distances(:) == inf), 68);	
			noninf =  distances(distances(:) ~= inf);
			obj.verifyEqual( sum(noninf),  2.826056638509035e+02, 'AbsTol', 1e-4);
		end
		
		%% Non symmetric connectivity.
		function non_symmetric_connectivity(obj)
			C = obj.linear_obj;
			
			C.set_connectivity_by_radius(3);
			conn = C.connectivity;
			conn(1:5:end,:);
			conn =  conn(conn(:,3) >= 0,:);
			
			C.connectivity = conn;
			
			% Finding a solution suffices
			C.shortest_path();
			
			C.curvature_penalty = 1e-100;
			C.shortest_path();
			
			C.torsion_penalty = 1e-100;
			C.shortest_path();
		end
		
		%% 
		function local_optimization(obj)
			import matlab.unittest.constraints.IsGreaterThanOrEqualTo;
			
			C = obj.linear_obj_large;
			C.length_penalty = 1;
			C.curvature_penalty = 1;

			[~,cost] = C.shortest_path();
			[~,lcost] = C.local_optimization();

			obj.verifyThat( cost.total, IsGreaterThanOrEqualTo(lcost.total) );
		end
		
		% Test that local optimization interploated new points.
		function local_optimization_interpolates_points(obj)
			C = obj.linear_obj;
			C.length_penalty = 1;
			C.curvature_penalty = 1;

			[curve1, cost] = C.shortest_path();
			num_points_1 = length(curve1);

			max_segment_length = 0.1;
			C.local_optimzation_max_curve_segment_length = max_segment_length;

			[curve2, lcost] = C.local_optimization();
			num_points_2 = length(curve2);
			obj.verifyGreaterThan(num_points_2, num_points_1);

			for i = 2:size(curve2, 1)
				d = norm(curve2(i-1, :) - curve2(i, :));
				obj.verifyLessThanOrEqual(d, 1.0001*max_segment_length);
			end

			C.local_optimzation_max_curve_segment_length = max_segment_length;
			[curve3, lcost] = C.local_optimization();
			num_points_3 = length(curve3);
			obj.verifyEqual(num_points_2, num_points_3);

			[curve4, cost] = C.shortest_path();
			num_points_4 = length(curve4);
			obj.verifyEqual(num_points_1, num_points_4);

			C.local_optimzation_max_curve_segment_length = max_segment_length;
			[curve5, lcost] = C.local_optimization();
			num_points_5 = length(curve5);
			obj.verifyEqual(num_points_2, num_points_5);
		end

		%% Constrained Shortest Paths
		function limits(obj)
			import matlab.unittest.constraints.*;
			C = obj.linear_obj_large;
			
			C.length_penalty = 0;
			C.torsion_penalty = 0;
			C.curvature_penalty = 0;
			
			length_max = 30;
			curvature_max = 5;		
			C.shortest_path();
	
			% Check that current solution does not fulfill limit
			% and enforce length.
			obj.verifyThat(C.info.length, IsGreaterThan(length_max));
			C.length_limit = length_max;
			C.shortest_path();		
			obj.verifyThat(C.info.length, IsLessThanOrEqualTo(length_max));
			
			% Enforce curvature.
			obj.verifyThat(C.info.curvature, IsGreaterThan(curvature_max));
			C.curvature_limit = curvature_max;
			C.shortest_path();
			obj.verifyThat(C.info.curvature, IsLessThanOrEqualTo(curvature_max));

			% Reset
			C.length_limit = inf;
			C.curvature_limit = inf;
			
			% Smaller problem for torsion
			C = obj.linear_obj;
			C.shortest_path();		
			torsion_max = 0.5;

			% Enforce torsion.
			obj.verifyThat(C.info.curvature, IsGreaterThan(torsion_max));
			C.torsion_limit = torsion_max;
			C.shortest_path();
			obj.verifyThat(C.info.curvature, IsLessThanOrEqualTo(torsion_max));
			
			C.torsion_limit = inf;
		end

		 %% Test the API for computing the visit tree.
 		function compute_visit_tree(obj)
			C = obj.linear_obj;
			C.set_connectivity_by_radius(3);
			C.length_penalty = 1.1;
			C.curvature_penalty = 0;
			C.torsion_penalty = 0;
			test_visit_tree(obj, C);
			C.curvature_penalty = 1e-100;
			test_visit_tree(obj, C);
			C.torsion_penalty = 1e-100;
			test_visit_tree(obj, C);
		end

		%
		function geodesic(obj)
			problem_size = [50 50];
			[xx,yy] = meshgrid(1:problem_size(1),1:problem_size(2));
			depth = sqrt((xx-25).^2 + (yy-25).^2);

			start_set = false(problem_size);
			end_set = false(problem_size);

			start_set(10,10) = true;
			end_set(end-10,end-10) = true;
			assert(depth(start_set) == depth(end_set));

			start_set(10,10) = true;
			end_set(end-10,end-10) = true;
			assert(depth(start_set) == depth(end_set));
			C = Curve_extraction('geodesic', depth, start_set, end_set);
			C.set_connectivity_by_radius(4);

			Cr = Curve_extraction('geodesic', depth, end_set,start_set);
			Cr.set_connectivity_by_radius(4);

			r = sqrt(2)*(25-10);
			max_dist = @(curve) max(sqrt((curve(:,1)-25).^2 + (curve(:,2)-25).^2));

			import matlab.unittest.constraints.*;

			for vd3 = [0.1 0.8 0.9 1 1.1 1.2 1.3 5 10];
				C.voxel_dimensions(3) = vd3;
				C.shortest_path();
				
				% Optimal curve for large enough height corresponds to great arc around the middle point.
				obj.verifyThat(max_dist(C.curve), IsLessThanOrEqualTo(r+1));
					
				% Start and stop set flipped.
				Cr.voxel_dimensions(3) = vd3;
				Cr.shortest_path();

				% Reverse problem have equal solution
				obj.verifyTrue(all(all(flipdim(Cr.curve,1) == C.curve)));
				
				C.local_optimization();
				obj.verifyThat(max_dist(C.curve), IsLessThanOrEqualTo(r+1));
			end
		end

		function explicit_data_term(obj)
			linear_data = rand(50,50);

			% Define local connectivity
			dims = ndims(linear_data);
			radius = 1; % 4-conn
			connectivity = get_all_directions(radius,dims);

			% Simple edge cost: 1/2 current + 1/2 target voxel value
			data = inf([size(linear_data) size(connectivity,1)]);
			problem_size = size(linear_data);

			for x = 1:size(data,1);
				for y = 1:size(data,2)
					for c= 1:size(connectivity,1);		
						tx = x + connectivity(c,1);
						ty = y + connectivity(c,2);
						
						if ((tx < 1) || (tx > problem_size(1)))
							continue;
						end
						
						if ((ty < 1) || (ty > problem_size(2)))
							continue;
						end
						
						data(x,y,c) = (linear_data(x,y)+ linear_data(tx,ty))/2;

					end
				end
			end

			%
			start_set = false(problem_size);
			end_set = false(problem_size);

			start_set(:,1) = true;
			end_set(:,end) = true;

			%% Compare to linear interpolation
			C = Curve_extraction('edge', data, connectivity, start_set, end_set);
			C.shortest_path();
			C1 = C.curve;

			%% The connectivity is chosen such that the edge cost will be exactly the same
			% as the explicitly given data costs.
			Cl = Curve_extraction('linear_interpolation', linear_data, start_set, end_set);
			Cl.set_connectivity_by_radius(radius);
			Cl.shortest_path();

			C2 = C.curve;
			obj.verifyTrue(all(C1(:) == C2(:)));
		end	
	end
end

