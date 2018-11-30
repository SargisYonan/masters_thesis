%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab

% https://www.mathworks.com/help/matlab/math/interpolating-scattered-data.html

clear all;
close all;
clc;

field_struct = load('generated_field.mat');
field = field_struct.field;
field_size = max(size(field));

figure();
top_pred_surf = pcolor(field); 
xlabel('$s_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$s_2$', 'Interpreter', 'latex', 'FontSize', 20);
shading interp; % gets rid of the grid lines on the surf()
title('Geospatially Autocorrelated Field', 'Interpreter', 'latex', 'FontSize', 20)
export_img_latex(gcf, '../figures/generated_field');

l=load('samples_save.mat');
X=l.X;
X = [X(1:7,:) ; X(9:end,:)]; % get rid of that weird one
sample_n = size(X,1);

hold on;

plot(X(:,1), X(:,2), 'ok', 'LineWidth', 2)

title('Sampled Geospatially Autocorrelated Field', 'Interpreter', 'latex', 'FontSize', 20)
export_img_latex(gcf, '../figures/sampled_generated_field');

%% Compute Variogram

Z = zeros(length(X),1);
for i = 1:length(X)
    Z(i) = field(X(i,1),X(i,2));
end

VAR = variogram([X(:,1), X(:,2)], Z, 'plot', false, 'type', 'gamma', 'nrbins', floor(length(Z)/2));
c0 = max(VAR.val);
a0 = max(VAR.distance)*2/3;

ix = linspace(0, max(VAR.distance));

kernels = {'gaussian'};

%%
figure()
plots = zeros(1, max(size(kernels)) + 1);

for i = 1:max(size(kernels))
    [range, sill, nugget, S] = variogramfit(VAR.distance, VAR.val, a0, c0, [], 'model', char(kernels(i)), 'plotit', false);
    b = [range, sill, nugget];
    kern = S.func(b, ix);
    plots(i) = plot(ix, kern, 'k-', 'LineWidth', 2);
    hold on
end

plot(VAR.distance, VAR.val, 'ro', 'LineWidth', 2);

title('Gaussian Variogram Fit To A Semi-Variogram', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('Lag $h$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('$\hat{\gamma}(h)$', 'Interpreter', 'latex', 'FontSize', 18)

axis([0, max(VAR.distance), 0, max(VAR.val)])


hold on
plot([0 range], [sill sill], ':k', 'LineWidth', 2)
hold on
plot([range range], [0 sill], '--r', 'LineWidth', 2)

legend('Gaussian Model', 'Empirical Semivariogram', 'sill', 'range', 'northwest', 'FontSize', 17)
grid on;

export_img_latex(gcf, '../figures/fit_kernel')

var_params = [range, sill, nugget];
variogram_model = S.func;

%% Calculate the covariance-variance matrix
n = length(Z);

C = zeros(n,n);
for i = 1:n
    x_t = X(i, 1);
    y_t = X(i, 2);
    for j = 1:n
        
        x_th = X(j, 1);
        y_th = X(j, 2);
        
        delta_x = x_t - x_th;
        delta_y = y_t - y_th;
        
        h = norm([delta_x delta_y]);
        
        if (h == 0)
            h = 0.1;
        end
        
        C(i,j) = variogram_model(var_params, h);
    end
end

pred_field = zeros(size(field));

%% Calculate approx. Cov() of each point
% calculate a vector, d, which contains
% all semiv values for an h = ||[i j] - [o_xi o_yi]||

% then compute a weighted-least squared to get that prediction value

d = zeros(n, 1);

warning('off','all')
warning

neighbor_conf = zeros(length(X), 2);

for i = 1:length(field)
    for j = 1:length(field)
        
        for k = 1:n
            x = X(k, 1);
            y = X(k, 2);

            delta_x = i - x;
            delta_y = j - y;

            h = norm([delta_x delta_y]);
            if (h == 0)
                h = 0.1;
            end
            
            d(k) = variogram_model(var_params, h);
        end
        
        pred_field(i,j) = Z' * (C\d);
        
%         idx = knnsearch(X, [i, j], 'k', 1);
%         neighbor_conf(idx, 1) = sum(d);
%         neighbor_conf(idx, 2) = neighbor_conf(idx, 2) + 1;
    end
end
% 
% for i = 1:length(neighbor_conf)
%     neighbor_conf(i, 1) = neighbor_conf(i, 1) / neighbor_conf(i, 2); % calculate average uncertainty of each neighbor
% end


figure()

pcolor(pred_field)
title('Kriging Predicted Field', 'Interpreter', 'latex', 'FontSize', 20)
shading interp; % gets rid of the grid lines on the surf()
xlabel('$s_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$s_2$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('\hat{Z}');

export_img_latex(gcf, '../figures/kriging_prediction');
% %%
% hold on
% 
% [vx,vy] = voronoi(X(:,1),X(:,2));
% 
% % Assign labels to the points.
% nump = length(X);
% plabels = arrayfun(@(n) {sprintf('N_{%d}\nv_{%d}=%.2f', n,n,100*neighbor_conf(n,1))}, (1:nump)');
% Hpl = text(X(:,1), X(:,2), plabels, 'FontWeight', ...
%       'bold', 'HorizontalAlignment','center');
% 
% %plot(X(:,1),X(:,2),'k+')
% %hold on;
% 
% title('Voronoi Tessellated Predicted Field', 'Interpreter', 'latex', 'FontSize', 20);
% plot(vx, vy, 'k-', 'LineWidth', 1.1)
% 
% export_img_latex(gcf, '../figures/natural_neighborhood_selection');
% 
% %%
% nn = 3;
% 
% edges = [];
% weights = [];
% for curr_n = 1:length(X)
% 
%     ix = knnsearch(X,X(curr_n,:),'k', nn + 1);
%     for j = 2:length(ix)
% 
%         if (length(edges) > 0)
%             if ( (sum(ismember(edges, [curr_n, ix(j)], 'rows')) == 0) && (sum(ismember(edges, [ix(j), curr_n], 'rows')) == 0) )
%                 edges = [edges; [curr_n, ix(j)]];
%                 weights = [weights ; 1/(abs(neighbor_conf(curr_n,1)) + abs(neighbor_conf(ix(j),1))) ];
%             end
%         else
%             edges = [edges; [curr_n, ix(j)]];
%             weights = [weights ; 1/(abs(neighbor_conf(curr_n,1)) + abs(neighbor_conf(ix(j),1))) ];
%         end
%     end
% end
% 
% vertex_names = arrayfun(@(n) {sprintf('N-%d', n)}, (1:length(X))');
% G = graph(edges(:,1)', edges(:,2)', weights', vertex_names);
% figure()
% 
% Gp = plot(G,'ko', 'LineWidth', 2)%, 'EdgeLabel', G.Edges.Weight)
% Gp.Marker = 'o';
% Gp.NodeColor = 'r';
% Gp.MarkerSize = 9;
% % turn off graph axes
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% title('Graph Representation of Voronoi Tessellated Field', 'Interpreter', 'latex', 'FontSize', 20)
% export_img_latex(gcf, '../figures/undir_graph')
% 
