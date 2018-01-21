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

%subplot(1,2,1);
top_pred_surf = pcolor(pred_field); 
shading interp; % gets rid of the grid lines on the surf()

%subplot(1,2,2);
%top_field_surf = pcolor(field); 
%shading interp; % gets rid of the grid lines on the surf()
        
%hold on;

%for s = 1:n
%    plot3(sampled_locations(s,1), sampled_locations(s,2), abs(samples(s)), ...
%        'ko', 'LineWidth',2,...
%        'MarkerEdgeColor','k',...
 %       'MarkerFaceColor',[0 0 0],...
  %      'MarkerSize',5)
  %  hold on;
%end
        
export_img_latex(gcf, 'idw_predicted_field');
save('idw_pred_field.mat', 'pred_field');
