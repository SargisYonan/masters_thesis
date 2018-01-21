%% Fun with Variogrpahy
clear all;
close all;
clc;

%% Load the test field and sample locations

X = load('test_field_10x10.mat');
X = X.field;

sample_locations = load('sample_locations_10x10.mat');
sample_locations = sample_locations.sample_locations;

field_size = max(size(X));
n = max(size(sample_locations));

observations = zeros(n, 1);

for i = 1:n
    x = sample_locations(i, 1);
    y = sample_locations(i, 2);
    observations(i) = X(x,y);
end

figure(1);
subplot(2,1,1);
surf(X)
title('Actual Field');

% create state vectors -- assumed to be passed in
states = [sample_locations, observations];

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

for i = 1:n^2
    % find all repeating instances of a given distance in the semivariogram    
    h = rep_dist(i);
    
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
 % plot the semivariogram
figure()
plot(distance, semivariogram, 'o')
hold on;

%% fit a 2nd degree polynomial to the semivarigram and plot it

ix = linspace(0,norm([field_size field_size]));
p = polyfit(distance, semivariogram, 2);
pfit = polyval(p, ix);
pfit_plot = plot(ix, pfit)
hold on;

%% calculate range and sill -- variogram model parameter
[sill, isx] = max(pfit);
range = ix(isx);

%% Plot the variogram

kern_plot = plot(ix, kernel(ix, range, sill, 'gaussian'));

grid on;
legend([pfit_plot, kern_plot], {'Least-Squares Fit Variogram', 'Kernel Model Fit Variogram'});

ylabel('Variogram');
xlabel('Distance');

%% Calculate the covariance-variance matrix
C = zeros(n,n);
for i = 1:n
    x_t = states(i, 1);
    y_t = states(i, 2);
    for j = 1:n
        x_th = states(j, 1);
        y_th = states(j, 2);
        
        delta_x = x_t - x_th;
        delta_y = y_t - y_th;
        
        h = norm([delta_x delta_y]);
        
        C(i,j) = kernel(h, range, sill, 'gaussian');
    end
end

%% For all possible positions, predict
% calculate a vector, d, which contains
% all semiv values for an h = ||[i j] - [o_xi o_yi]||

% then compute a weighted-least squared to get that prediction value

pred_X = zeros(field_size,field_size);
d = zeros(n, 1);
for i = 1:field_size
    for j = 1:field_size
        
        for k = 1:n
            x = states(k, 1);
            y = states(k, 2);

            delta_x = i - x;
            delta_y = j - y;

            h = norm([delta_x delta_y]);

            d(k) = kernel(h, range, sill, 'gaussian');
        end
        
        pred_X(i,j) = observations' * (C\d);
        
    end
end

figure(1);
subplot(2,1,2);
surf(pred_X)
title('Predicted Field')

