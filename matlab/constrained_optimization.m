
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

hold on;
xlim([sigma1 sigma2]);
xlabel('\sigma')
ylabel('Dual function value.')
drawnow;

c = [-1; 0];
% A = [1 -g1;
%      1 -g2];
% b = [d1 - g1*sigma1;
%      d2 - g2*sigma2];
A = zeros(0, 2);
b = zeros(0, 1);
lb = [-inf; sigma1];
ub = [inf;  sigma2];

for iter = 1:20
	if iter == 1
		sigma = mean([sigma1 sigma2]);
		projected_energy = nan;
	else
		[x, fval, exitflag] = linprog(c, A, b, [], [], lb, ub);
		assert(exitflag == 1);

		projected_energy = x(1);
		sigma            = x(2);
		plot(sigma, projected_energy, 'k*')
		if iter == 2
			ylim([0 1.1*projected_energy]);
		end
	end

	C.curvature_penalty = sigma;
	[Curve, Enew] = C.shortest_path();
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
plot(Curve(:,2), Curve(:,1),'color', 'red', 'linewidth',5)
plot(Curve(:,2), Curve(:,1),'color', 'red', 'linewidth',5)
plot(Curve(:,2), Curve(:,1),'color', 'red', 'linewidth',5)
title('Length regularization with curvature constraint.');


