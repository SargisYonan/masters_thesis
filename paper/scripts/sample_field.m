%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

plot_sample_points = true;

field_struct = load('generated_field.mat');
field = field_struct.field;
field_size = max(size(field));
number_of_samples = 50;

sample_locations = zeros(number_of_samples, 2);
samples = zeros(number_of_samples, 1);

% initialize random
seed = floor(mod(now * 1E6, 1000));
rng(seed, 'twister');

if plot_sample_points
    top_surf = pcolor(field); 
    xlabel('x_1');
    ylabel('x_2');
    zlabel('u');
    shading interp; % gets rid of the grid lines on the surf()
    hold on;
end

for i = 1:number_of_samples
    sample_locations(i, :) = randi([1 field_size], 1, 2);
    samples(i) = field(sample_locations(i,1), sample_locations(i,2));
    
    if plot_sample_points
        plot3(sample_locations(i,1), sample_locations(i,2), abs(samples(i)), ...
            'ko', 'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 0],...
            'MarkerSize',5)
        hold on;
    end
    
end

xlabel('x_1');
ylabel('x_2');
zlabel('u');

if plot_sample_points
   export_img_latex(gcf, '../figures/sampled_generated_field')
end

save('sampled_points.mat', 'samples');
save('sampled_locations.mat', 'sample_locations');




