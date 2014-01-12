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

length_reg =  1;

LAB = @(I) double(applycform(uint8(I), makecform('srgb2lab')));
RGB = @(I) double(applycform(uint8(I),makecform('lab2srgb')));

% Image
I = imread('../data/irrawaday-delta-sevcik.jpg');
I = I(40:end,:,:);
I = imresize(I,0.25);

% mesh_map defining allowed pixels
% Encoded as
% 0: Disallowed
% 1: Allowed
% 2: Start set
% 3: End set.
[w,h,c] = size(I);
mesh_map = ones([size(I,1) size(I,2)], 'int32');

% Define start/end set.
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

unary = inf(size(I,1), size(I,2));
target_intensity = zeros(size(I));
for i = 1 : numel(River_colors);
    
    target_intensity(:,:,1) = River_colors{i}(1);  % R
    target_intensity(:,:,2) = River_colors{i}(2);  % G
    target_intensity(:,:,3) = River_colors{i}(3);  % B
    
    unary_i = LAB(I) - LAB(target_intensity);
    unary_i = sum(unary_i.^2, 3);
    unary = min(unary, unary_i);
end

% Scale down
unary = unary/1e3;

% Setup the problem instance
C = Curve_extraction(unary, start_set, end_set);

% Use all edges with size <= regularization radius
C.set_connectivity_by_radius(2.5);


%% Define settings
% Use all edges with size <= regularization radius
C.regularization_radius = 2.5;
C.verbose = true;

% rho in paper
C.length_penalty = 1;

% sigma in paper
C.curvature_penalty = 0;

% tau in paper
C.torsion_penalty = 0;
C.use_a_star = false;

%% Solve

tree = C.compute_tree()

% Placeholder "visualization."
imagesc(tree);




