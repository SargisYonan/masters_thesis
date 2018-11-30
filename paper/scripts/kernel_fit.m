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
plots(max(size(kernels)) + 1) = plot(VAR.distance, VAR.val, 'ko', 'LineWidth', 1);
hold on;

i = 1;
[range, sill, nugget, S] = variogramfit(VAR.distance, VAR.val, a0, c0, [], 'model', char(kernels(i)), 'plotit', false);
b = [range, sill, nugget];
kern = S.func(b, ix);
plots(i) = plot(ix, kern, ':r', 'LineWidth', 3);
hold on

i = 2;
[range, sill, nugget, S] = variogramfit(VAR.distance, VAR.val, a0, c0, [], 'model', char(kernels(i)), 'plotit', false);
b = [range, sill, nugget];
kern = S.func(b, ix);
plots(i) = plot(ix, kern, '--k', 'LineWidth', 3);
hold on

i = 3;
[range, sill, nugget, S] = variogramfit(VAR.distance, VAR.val, a0, c0, [], 'model', char(kernels(i)), 'plotit', false);
b = [range, sill, nugget];
kern = S.func(b, ix);
plots(i) = plot(ix, kern, '-b', 'LineWidth', 3);
hold on

title('Variogram Model Comparison', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('Lag $h$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('$\gamma(h)$', 'Interpreter', 'latex', 'FontSize', 18)

axis([0, max(VAR.distance), 0, max(VAR.val)])
legend(plots, 'Gaussian Model', 'Spherical Model', 'Exponential Model', 'Experimental Semivariogram', 'Location', 'northwest')

export_img_latex(gcf, '../figures/fit_kern_comp')
