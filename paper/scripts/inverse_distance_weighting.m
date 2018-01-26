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


pred_field = zeros(field_size, field_size);
p = 3;

for i = 1:field_size
    for j = 1:field_size
        
        pred = 0.0;
        w = 0.0;
        w_sum = 0.0;
        
        for k = 1:n
            h = norm([i j] - sampled_locations(k,:));
            
            if (h == 0) % just use the observed value
                pred = samples(k);
                w_sum = 1.0;
            else
                w = 1/(h^p);
                w_sum = w_sum + w;
                pred = pred + w*(samples(k));
            end
            
        end
        pred = pred / w_sum;
        
        pred_field(i,j) = pred;
    end
end

figure(1);
top_pred_surf = pcolor(pred_field); 
xlabel('x_1');
ylabel('x_2');
zlabel('u');
shading interp; % gets rid of the grid lines on the surf()

export_img_latex(gcf, 'idw_predicted_field');
save('idw_pred_field.mat', 'pred_field');

figure(2);
surf(pred_field); 
xlabel('x_1');
ylabel('x_2');
zlabel('u');
export_img_latex(gcf, 'idw_side_pred_field');

