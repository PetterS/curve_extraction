% This script reproduces Figure 6 in
% Shortest Paths with Curvature and Torsion
% Petter Strandmark, Johannes Ulén, Fredrik Kahl, Leo Grady. 
% International Conference on Computer Vision. 2013.
%
% Windows note:
% If you are using precompiled mex files
% you need Visual C++ Redistributable Packages for Visual Studio 2013:
% http://www.microsoft.com/en-us/download/details.aspx?id=40784
close all
clear all
addpath([fileparts(mfilename('fullpath')) filesep '..']);

% Regularizations
use_a_star = true;
curvature_regs = [0.25 50 500];
length_regs = [0.25 1 1000];

LAB = @(I) double(applycform(uint8(I), makecform('srgb2lab')));
RGB = @(I) double(applycform(uint8(I),makecform('lab2srgb')));

% Image
I = imread('../../data/irrawaday-delta-sevcik.jpg');
I = I(40:end,:,:);
I = imresize(I,0.25);

% Define start/end set.
[w,h,c] = size(I);
start_set = false([size(I,1) size(I,2)]);
end_set = false([size(I,1) size(I,2)]);

start_set(1:3, 1:round(end/2)) = 2;
start_set(end,1:round((1/20*h))) = 2;
end_set(end-3:end, :) = 3;

% Remove overlap
start_set(end_set) = 0;

%% Unary, hand crafted
River_colors = {};
River_colors{end+1} = [169 151 137];
River_colors{end+1} = [200 156 119];
River_colors{end+1} = [255 255 255];
River_colors{end+1} = [191 194 202];
River_colors{end+1} = [233 225 233];
River_colors{end+1} = [131 116 83];
River_colors{end+1} = [182 149 121];
River_colors{end+1} = [205 176 165];
River_colors{end+1} = [156 133 102];

data = inf(size(I,1), size(I,2));
target_intensity = zeros(size(I));
for i = 1 : numel(River_colors);
    
    target_intensity(:,:,1) = River_colors{i}(1);  % R
    target_intensity(:,:,2) = River_colors{i}(2);  % G
    target_intensity(:,:,3) = River_colors{i}(3);  % B
    
    data_i = LAB(I) - LAB(target_intensity);
    data_i = sum(data_i.^2, 3);
    data = min(data, data_i);
end

% Scale down
data = data/1e3;

% Setup the problem instance
C = Curve_extraction(data, start_set, end_set);

C.num_threads = min(feature('numThreads'),2);
C.descent_method = 'lbfgs';

% Use all edges with size <= regularization radius
C.set_connectivity_by_radius(2.5);


%% Define settings
% Use all edges with size <= regularization radius
C.verbose = true;

% rho in paper
C.length_penalty = 0;

% sigma in paper
C.curvature_penalty = 0;

% tau in paper
C.torsion_penalty = 0;
C.use_a_star = use_a_star;

%% 
m = [0.0384 0.1112 0.2790;
    0.6650 0.1675 0.1675;
    0.3726 0.9210 0.1350];

% % Trying different lengths
assert(length(curvature_regs) == 3);

C.length_penalty = 0;
C.curvature_penalty = curvature_regs(1);
curvature_1 = C.shortest_path();

C.curvature_penalty = curvature_regs(2);
curvature_2 = C.shortest_path();

C.curvature_penalty = curvature_regs(3);
curvature_3 = C.shortest_path();

% Display


%%
assert(length(length_regs) == 3);

C.curvature_penalty = 0;
C.length_penalty = length_regs(1);
length_1 = C.shortest_path();

C.length_penalty = length_regs(2);
length_2 = C.shortest_path();

C.length_penalty = length_regs(3);
length_3 = C.shortest_path();

%% Generate Figure 6 (c) and (d) visit order of the algorithm
% with and without A*.
C.curvature_penalty = curvature_regs(2);
C.length_penalty = 0;
C.store_visit_time = true;

C.use_a_star = false;
C.shortest_path();
without_astar = double(C.visit_map);

C.use_a_star = true;
C.shortest_path();
with_astar = double(C.visit_map);

%% Display
figure(1);
clf;imagesc(rgb2gray(I)); 
colormap gray(256); hold on; axis equal; hold on; axis off;
plot(length_1(:,2), length_1(:,1),'color', m(1,:) , 'linewidth',5)
plot(length_2(:,2), length_2(:,1),'color', m(2,:) , 'linewidth',5)
plot(length_3(:,2), length_3(:,1),'color', m(3,:) , 'linewidth',5)
title('Figure 6 (a): Length regularization.');

%%
figure(2);
clf;
imagesc(rgb2gray(I))
colormap gray(256); hold on; axis equal; axis off;
plot(curvature_1(:,2),curvature_1(:,1),'color', m(1,:) , 'linewidth',5)
plot(curvature_2(:,2),curvature_2(:,1),'color', m(2,:) , 'linewidth',5)
plot(curvature_3(:,2),curvature_3(:,1),'color', m(3,:) , 'linewidth',5)
title('Figure 6 (b): Curvature regularization.');

%%
figure(3); clf; 
Ibw = double(rgb2gray(I));

unvisited = (without_astar == -1);

without_astar_m = without_astar+max(Ibw(:));
without_astar_m(unvisited) = Ibw(unvisited);
imagesc(without_astar_m);

colormap([gray(max(Ibw(:))); jet(max(without_astar_m(:)))]);

% Windows cannot render colormaps with more than 256 colors without OpenGL.
if (ispc)
	set(gcf,'renderer','OpenGL')
end

title('Figure 6 (c): Visit order for medium curvature regularization \bf{without A*}. ')
axis equal; axis off;

%%
figure(4); clf;
unvisited = (with_astar == -1);

with_astar_m = with_astar+max(Ibw(:));
with_astar_m(unvisited) = Ibw(unvisited);
imagesc(with_astar_m);

colormap([gray(max(Ibw(:))); jet(max(with_astar_m(:)))]);
if (ispc)
	set(gcf,'renderer','OpenGL')
end

title('Figure 6 (d):  Visit order for medium curvature regularization \bf{with A*}.')
axis equal; axis off;


