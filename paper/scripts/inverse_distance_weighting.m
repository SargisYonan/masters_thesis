%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

field_struct = load('generated_field.mat');
field = field_struct.field;
field_size = max(size(field));
% 
% sample_struct = load('sampled_points.mat');
% samples = sample_struct.samples;
% n = max(size(samples));

l=load('samples_save.mat');
sampled_locations = l.X;
sampled_locations = [sampled_locations(1:7,:) ; sampled_locations(9:end,:)]; % get rid of that weird one
n = size(sampled_locations,1);

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
                pred = field(sampled_locations(k,1),sampled_locations(k,2));
                w_sum = 1.0;
            else
                w = 1/(h^p);
                w_sum = w_sum + w;
                pred = pred + w * field(sampled_locations(k,1),sampled_locations(k,2));
            end
            
        end
        pred = pred / w_sum;
        
        pred_field(i,j) = pred;
    end
end

figure(1);
top_pred_surf = pcolor(pred_field); 
shading interp; % gets rid of the grid lines on the surf()

xlabel('$s_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$s_2$', 'Interpreter', 'latex', 'FontSize', 20);

title('Inverse Distance Weighting Prediction', 'Interpreter', 'latex', 'FontSize', 20)
export_img_latex(gcf, '../figures/idw_predicted_field');
save('idw_pred_field.mat', 'pred_field');

figure(2);
surf(pred_field); 
xlabel('$s_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$s_2$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$\hat{Z}$', 'Interpreter', 'latex', 'FontSize', 20);
export_img_latex(gcf, 'idw_side_pred_field');
