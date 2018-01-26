%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

kernel_model = 'spherical';

F = load('generated_field.mat');
field = F.field;
field_size = max(size(field));

S = load('sampled_points.mat');
observations = S.samples;
n = max(size(observations));

O = load('sampled_locations.mat');
observation_locations = O.sample_locations;

D = load('kernel_params.mat');
semivariogram = D.semivariogram;
distances = D.distance;



VAR = variogram([observation_locations(:,1), observation_locations(:,2)], observations, 'plot', false, 'type', 'gamma', 'nrbins', floor(max(size(observations))/2));
c0 = max(VAR.val);
a0 = max(VAR.distance)*2/3;
[range, sill, nugget, VF] = variogramfit(VAR.distance, VAR.val, a0, c0, 1, 'model', kernel_model, 'solver','fminsearchbnd', 'plotit', false);
var_params = [range, sill, nugget];
variogram_model = VF.func;

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
        
        if (h == 0)
            h = 0.1;
        end
        
        C(i,j) = variogram_model(var_params, h);%kernel(h, range, sill, kernel_model);
    end
end

pred_field = zeros(field_size, field_size);

%% Calculate approx. Cov() of each point
% calculate a vector, d, which contains
% all semiv values for an h = ||[i j] - [o_xi o_yi]||

% then compute a weighted-least squared to get that prediction value

d = zeros(n, 1);

warning('off','all')
warning

for i = 1:field_size
    for j = 1:field_size
        
        for k = 1:n
            x = observation_locations(k, 1);
            y = observation_locations(k, 2);

            delta_x = i - x;
            delta_y = j - y;

            h = norm([delta_x delta_y]);
            if (h == 0)
                h = 0.1;
            end
            
            d(k) = variogram_model(var_params, h); %kernel(h, range, sill, kernel_model);
        end
        
        pred_field(i,j) = observations' * (C\d);
        
    end
end

figure(1);
surf(pred_field)
title('Predicted Field')
xlabel('x_1');
ylabel('x_2');
zlabel('u');
export_img_latex(gcf, 'kriging_side_pred_field');

figure(2);
surf(field)
title('Actual Field')
xlabel('x_1');
ylabel('x_2');
zlabel('u');

figure(3);
pcolor(pred_field)
title('Predicted Field')
shading interp; % gets rid of the grid lines on the surf()
xlabel('x_1');
ylabel('x_2');
zlabel('u');
export_img_latex(gcf, 'kriging_top_pred_field');

figure(4);
pcolor(field)
title('Actual Field')
shading interp; % gets rid of the grid lines on the surf()
xlabel('x_1');
ylabel('x_2');
zlabel('u');

