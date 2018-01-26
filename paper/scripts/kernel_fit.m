%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

F = load('generated_field.mat');
field = F.field;
field_size = max(size(field));

S = load('sampled_points.mat');
observations = S.samples;
n = max(size(observations));

O = load('sampled_locations.mat');
observation_locations = O.sample_locations;

VAR = variogram([observation_locations(:,1), observation_locations(:,2)], observations, 'plot', false, 'type', 'gamma', 'nrbins', floor(max(size(observations))/2));
c0 = max(VAR.val);
a0 = max(VAR.distance)*2/3;

ix = linspace(0, max(VAR.distance));

kernels = {'gaussian', 'spherical', 'exponential'};

figure()
plots = zeros(1, max(size(kernels)) + 1);

for i = 1:max(size(kernels))
    [range, sill, nugget, S] = variogramfit(VAR.distance, VAR.val, a0, c0, [], 'model', char(kernels(i)), 'plotit', false);
    b = [range, sill, nugget];
    kern = S.func(b, ix);
    plots(i) = plot(ix, kern);
    hold on
end

plots(max(size(kernels)) + 1) = plot(VAR.distance, VAR.val, 'k.');

title('Variogram')
xlabel('Distance')
ylabel('\gamma')

axis([0, max(VAR.distance), 0, max(VAR.val)])
legend(plots, 'Gaussian Model', 'Spherical Model', 'Exponential Model', 'Experimental Variogram', 'Location', 'northwest')

export_img_latex(gcf, 'fit_kern_comp.png')
