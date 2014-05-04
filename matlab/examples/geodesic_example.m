close all; clear all;
addpath([fileparts(mfilename('fullpath')) filesep '..']);

problem_size = [50 50];
[xx,yy] = meshgrid(1:problem_size(1),1:problem_size(2));
depth = 10 * sqrt((xx-25).^2 + (yy-25).^2);

start_set = false(problem_size);
end_set = false(problem_size);

start_set(10,10) = true;
end_set(end-10,end-10) = true;
assert(depth(start_set) == depth(end_set));

r = sqrt(2)*(25-10);
m = [25 25];
x = r*cos(linspace(0,2*pi,100))+m(1);
y = r*sin(linspace(0,2*pi,100))+m(2);

C = Geodesic_shortest_path(depth, start_set, end_set);
C.set_connectivity_by_radius(5);


%% Shortest path
figure(1);
C.shortest_path();
C.plot_curve();
title(sprintf('Curve length: %g\nHalf-circle length: %g\n', C.length(), r*pi));
plot(x,y,'-.');

%% Local optimization
figure(2);
C.local_optimzation_max_curve_segment_length = 0.5;
C.local_optimization;
C.plot_curve;
title(sprintf('Curve length: %g\nHalf-circle length: %g\n', C.length(), r*pi));
plot(x,y,'-.');




