% Showcasing the shortest_path tree function;
close all
clear all

rng(1)
data = double(imread('cameraman.tif'));

% Quiver plot is extremly demanding even on a modern computers
% solves a for a smaller image.
data = imresize(data,0.25);


% True at voxels where it allowed to start/end
start_set = false(size(data));
end_set = false(size(data));

start_set(64,11) = true;
end_set(64,52) = true;

%% Create Curve object
data_type = 'linear_interpolation';
C = Curve_extraction(data_type, data, start_set, end_set);
C.set_connectivity_by_radius(4);



%% Length
C.length_penalty = 1;
C.curvature_penalty = 0;
[tree,curve] = C.compute_visit_tree();
C.curve = curve;

figure(1);
C.display();
hold on;

X = [];
Y = [];
U = [];
V = [];
for i = 1:size(tree, 1)
    for j = 1:size(tree, 2)
        if tree(i, j) > 0
            [i2, j2] = ind2sub(size(tree), tree(i, j));
            X(end+1) = i;
            Y(end+1) = j;
            U(end+1) = i2 - i;
            V(end+1) = j2 - j;
        end
    end
end
quiver(Y,X,U,V, 0);
axis ij; axis equal;
title('Shortest-path tree length regularization');


%% Curvature
figure(2);
C.length_penalty = 0;
C.curvature_penalty = 1;
[tree,curve] = C.compute_visit_tree();
C.curve = curve;
C.display();
hold on;

X = [];
Y = [];
U = [];
V = [];
for i = 1:size(tree, 1)
    for j = 1:size(tree, 2)
        if tree(i, j) > 0
            [i2, j2] = ind2sub(size(tree), tree(i, j));
            X(end+1) = i;
            Y(end+1) = j;
            U(end+1) = i2 - i;
            V(end+1) = j2 - j;
        end
    end
end
quiver(Y,X,U,V, 0);
axis ij; axis equal;
title('Shortest-path tree curvature regularization');