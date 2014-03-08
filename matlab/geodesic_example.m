clc; close all; clear all;
addpath([fileparts(mfilename('fullpath')) filesep 'library']);

problem_size = [50 50];
[xx,yy] = meshgrid(1:problem_size(1),1:problem_size(2));
depth = sqrt((xx-25).^2 + (yy-25).^2);

start_set = false(problem_size);
end_set = false(problem_size);

start_set(10,10) = true;
end_set(end-10,end-10) = true;
assert(depth(start_set) == depth(end_set));

r = sqrt(2)*(25-10);
m = [25 25];
x = r*cos(linspace(0,2*pi,100))+m(1);
y = r*sin(linspace(0,2*pi,100))+m(2);

C = Curve_extraction('geodesic', depth, start_set, end_set);

% voxel_dimenions(3):
% A small number leads to a almost flat surface
% A very large number gives a shortest pathfollowing the 
% level sets of the surface.
C.voxel_dimensions(3) = 100;
C.set_connectivity_by_radius(4);
C.length_penalty = 5; 
C.curvature_penalty = 0;

%% Shortest path
figure(1);
C.shortest_path();
C.plot_curve();
plot(x,y,'-.');

%% Local optimization
figure(2);
C.local_optimization;
C.plot_curve;
plot(x,y,'-.');