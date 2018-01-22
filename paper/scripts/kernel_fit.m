%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

clear all;
close all;
clc;

D = load('kernel_params.mat');
range = D.range;
sill = D.sill;
semivariogram = D.semivariogram;
distances = D.distance;

ix = linspace(0, max(distances));

figure()

gauss_kern = kernel(ix, range/2, sill, 'gaussian');
plot(ix, gauss_kern, 'r')
hold on

sph_kern = kernel(ix, range/2, sill, 'sph');
plot(ix, sph_kern, 'k')
hold on

legend('Gaussian Kernel', 'Spherical Kernel')

plot(distances, semivariogram, 'bo')

export_img_latex(gcf, 'fit_kern_comp.png')
