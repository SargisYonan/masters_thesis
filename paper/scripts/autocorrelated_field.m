%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

%% Create a random gaussian distributed autocorrelated field
% and save to a file
% field settings
field_size = 100; % final area will be this square area
field_distribution = 'gaussian';
field_gain = 100;

[X,~] = meshgrid(0: field_size - 1);
Z = randn(size(X)); % normally distributed -- used rand() for a single round point

sigma = 8;%32; % inverse noise
% autocorrelation
field = field_gain * imfilter(Z, fspecial(field_distribution, [field_size, field_size], sigma));

save('diff_generated_field.mat', 'field');

%% create a plot for the paper

figure(1);
top_surf = pcolor(field); 
xlabel('x_1');
ylabel('x_2');
zlabel('u');
shading interp; % gets rid of the grid lines on the surf()

% export_img_latex(gcf, 'generated_field_top_view')


figure(2);
side_surf = surf(field);
xlabel('x_1');
ylabel('x_2');
zlabel('u');

% export_img_latex(gcf, 'generated_field_side_view')
