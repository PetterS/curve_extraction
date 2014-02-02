
close all
clear all
clc

% Curvature must be less than this.
curvature_limit = 1;

% Lower and upper bounds on regularization.
sigma1 = 0.25;
sigma2 = 1000;

LAB = @(I) double(applycform(uint8(I), makecform('srgb2lab')));
RGB = @(I) double(applycform(uint8(I),makecform('lab2srgb')));

% Image
I = imread('../data/irrawaday-delta-sevcik.jpg');
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
data_type = 'linear_interpolation';
C = Curve_extraction(data_type, data, start_set, end_set);

C.num_threads = min(feature('numThreads'),2);
C.descent_method = 'lbfgs';

% Use all edges with size <= regularization radius
C.set_connectivity_by_radius(2.5);


%% Define settings
% Use all edges with size <= regularization radius
C.verbose = true;

% rho in paper
C.length_penalty = 0.05;

% sigma in paper
C.curvature_penalty = 0;

% tau in paper
C.torsion_penalty = 0;
C.use_a_star = false;

%% 
m = [0.0384 0.1112 0.2790;
    0.6650 0.1675 0.1675;
    0.3726 0.9210 0.1350];



hold on;
xlim([sigma1 sigma2]);
xlabel('\sigma')
ylabel('Dual function value.')
drawnow;

C.curvature_penalty = sigma1;
[~, E1] = C.solve();
C1 = E1.curvature / sigma1; % Integrated curvature.
g1 = C1 - curvature_limit;  % Supergradient.
d1 = E1.total - sigma1*curvature_limit; % Dual function value.

C.curvature_penalty = sigma2;
[~, E2] = C.solve();
C2 = E2.curvature / sigma2; % Integrated curvature.
g2 = C2 - curvature_limit;  % Supergradient.
d2 = E2.total - sigma2*curvature_limit; % Dual function value.

plot(sigma1, d1, 'k.');
plot(sigma2, d2, 'k.');
draw_line(sigma1, d1, g1, sigma1, sigma2, 'k-');
draw_line(sigma2, d2, g2, sigma1, sigma2, 'k-');
drawnow;

%% 
assert(g1 > 0 && g2 < 0);

c = [-1; 0];
A = [1 -g1;
     1 -g2];
b = [d1 - g1*sigma1;
     d2 - g2*sigma2];
gap = 0;

for iter = 1:20
	[x, fval, exitflag] = linprog(c, A, b);
	assert(exitflag == 1);

	projected_energy = x(1);
	sigma            = x(2);
	plot(sigma, projected_energy, 'k*')
	if iter == 1
		ylim([0 1.1*projected_energy]);
	end

	C.curvature_penalty = sigma;
	[Curve, Enew] = C.solve();
	Cnew = Enew.curvature / sigma;
	gnew = Cnew - curvature_limit;
	dnew = Enew.total - sigma*curvature_limit;
	
	gap = projected_energy - dnew;
	
	A = [A; 1, -gnew];
	b = [b; dnew - gnew*sigma];
	
	plot(sigma, dnew, 'k.');
	plot([sigma, sigma], [projected_energy, dnew], 'k:');
	draw_line(sigma, dnew, gnew, sigma1, sigma2, 'k-');
	drawnow;
	
	fprintf('sigma=%g, gap=%g\n', sigma, gap);
	if gap < 1e-3
		break;
	end
end

fprintf('Final curvature: %g\n', Cnew);

figure;
imagesc(rgb2gray(I)); 
colormap gray(256); hold on; axis equal; hold on; axis off;
plot(Curve(:,2), Curve(:,1),'color', m(1,:) , 'linewidth',5)
plot(Curve(:,2), Curve(:,1),'color', m(2,:) , 'linewidth',5)
plot(Curve(:,2), Curve(:,1),'color', m(3,:) , 'linewidth',5)
title('Length regularization with curvature constraint.');


