%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

kernel_model = 'sph';

K = load('kernel_params.mat');
range = K.range;
sill = K.sill;

F = load('generated_field.mat');
field = F.field;
field_size = max(size(field));

S = load('sampled_points.mat');
observations = S.samples;
n = max(size(observations));

O = load('sampled_locations.mat');
observation_locations = O.sample_locations;

%% Calculate the covariance-variance matrix
C = zeros(n,n);
for i = 1:n
    x_t = observation_locations(i, 1);
    y_t = observation_locations(i, 2);
    for j = 1:n
        
        x_th = observation_locations(j, 1);
        y_th = observation_locations(j, 2);
        
        delta_x = x_t - x_th;
        delta_y = y_t - y_th;
        
        h = norm([delta_x delta_y]);
        
        C(i,j) = kernel(h, range, sill, kernel_model);
    end
end

pred_field = zeros(field_size, field_size);

%% Calculate approx. Cov() of each point
% calculate a vector, d, which contains
% all semiv values for an h = ||[i j] - [o_xi o_yi]||

% then compute a weighted-least squared to get that prediction value

d = zeros(n, 1);

for i = 1:field_size
    for j = 1:field_size
        
        for k = 1:n
            x = observation_locations(k, 1);
            y = observation_locations(k, 2);

            delta_x = i - x;
            delta_y = j - y;

            h = norm([delta_x delta_y]);

            d(k) = kernel(h, range, sill, kernel_model);
        end
        
        pred_field(i,j) = observations' * (C\d);
        
    end
end

figure(1);
subplot(2,1,1);
surf(pred_field)
title('Predicted Field')
subplot(2,1,2);
surf(field)
title('Actual Field')

figure(2);
subplot(2,1,1);
pcolor(pred_field)
title('Predicted Field')
subplot(2,1,2);
pcolor(field)
title('Actual Field')

