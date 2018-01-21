%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

field_struct = load('generated_field.mat');
field = field_struct.field;
field_size = max(size(field));

sample_struct = load('sampled_points.mat');
samples = sample_struct.samples;
n = max(size(samples));

sample_loc_struct = load('sampled_locations.mat');
sampled_locations = sample_loc_struct.sample_locations;

% create state vectors -- assumed to be passed in
states = [sampled_locations, samples];

%% calculate state respective means

mu_x = mean(states(:,1));
mu_y = mean(states(:,2));
mu_z = mean(states(:,3));
sigma_z = std(states(:,3));

mu = [mu_x mu_y mu_z];

%% create covariance-variance matrix

% constructing a semi-variogram, using the fact that
% E(Z_t - Z_t+h)^2 = 2\semigamma(h) is a variogram
% NOTE: E(Z_t - Z_t+h)^2 = normpdf((Z_t - Z_t+h), mu_z, sigma_z)^2
% -- h is the lag distance = ||X-Y||

rep_dist = [];
sem = [];

for i = 1:n
    x_t = states(i, 1);
    y_t = states(i, 2);
    z_t = states(i, 3);
    
    for j = 1:n
        x_th = states(j, 1);
        y_th = states(j, 2);
        z_th = states(j, 3);
        
        delta_x = x_t - x_th;
        delta_y = y_t - y_th;
        
        h = norm([delta_x delta_y]);
        
        rep_dist(end + 1) = h;
        sem(end + 1) = ( (z_t - z_th) * normpdf((z_t - z_th), mu_z, sigma_z) )^2;
    end
end

%% compute the semivariogram
% eq. 3 on http://www.stat.ucla.edu/~nchristo/statistics_c173_c273/cressie_85a.pdf
% where 2\gamma(h) = Sum_1-Nh(\semigamma(h))

semivariogram = [];
distance = [];
max_var_dist = norm([field_size field_size])/4;

for i = 1:n^2
    % find all repeating instances of a given distance in the semivariogram    
    h = rep_dist(i);
    if (h < max_var_dist)
        % you only need to average the semivariogram for the given distance, h,
        % if you havent previously
        if (~any(distance == h))
            ix = find(rep_dist == h);

            % average the values in the semivariogram and construct 
            % final semivariogram and distance arrays
            semivariogram(end + 1) = mean(sem(ix));
            distance(end + 1) = h;
        end
    end
end

% plot the semivariogram 
figure()
semi_v_plot = plot(distance, semivariogram, ...
    'ko', 'LineWidth', 2,...
    'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', [0 0 0],...
    'MarkerSize', 2);
hold on;

%% fit a 2nd degree polynomial to the semivarigram and plot it

ix = linspace(1, norm([field_size field_size]));
p = polyfit(distance, semivariogram, 2);
pfit = polyval(p, ix);
pfit_plot = plot(ix, pfit);
axis([0 max_var_dist 0 max(semivariogram)])
hold on;

ylabel('Variogram');
xlabel('Distance');
legend([pfit_plot, semi_v_plot], {'Least-Squares Fit Variogram', 'Experimental Variogram'});

%% calculate range and sill -- variogram model parameter

[sill, isx] = max(pfit);
range = ix(isx);

%% pack it up
save('kernel_params.mat', 'sill', 'range');
