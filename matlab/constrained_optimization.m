close all; clear all;
rng(0)
n = 100;
data = 100*rand(n,n);

start_set = false(size(data));
end_set = false(size(data));

% Define start/end set.
start_set(1:25,1) = true;
end_set(end-25:end,end) = true;

data_type = 'linear_interpolation';
C = Curve_extraction(data_type, data, start_set, end_set);
C.set_connectivity_by_radius(3);

C.length_penalty = 0.25;
C.curvature_penalty = 0.25;

figure(1);

C.shortest_path();
C.plot_curve()

figure(2);
% Upper bounds on the length and curvature.
C.length_limit = 115;
C.curvature_limit = 2;

C.shortest_path();
C.plot_curve()