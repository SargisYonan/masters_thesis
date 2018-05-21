%% SARGIS S YONAN
% Master Thesis
% Autonomous Systems Lab
% Uncertaianty Suppression
% https://www.mathworks.com/help/matlab/math/interpolating-scattered-data.html

clear all;
close all;
clc;

field_struct = load('new_generated_field.mat');
field = field_struct.field;
field_size = max(size(field));

stream = RandStream.getGlobalStream;
RandStream.seed = now;

n = 5;

% take a random measurments from [1,n]
X = zeros(n, 2);
X(:,1) = randi(length(field), n, 1);
X(:,2) = randi(length(field), n, 1);

max_runs = 20;
uncerts = zeros(max_runs,1);
errors = zeros(max_runs, 1);
sample_track = zeros(max_runs, 1);

for runs = 1:max_runs

%% Compute Variogram
n = size(X,1);
Z = zeros(n, 1);
for i = 1:n
    Z(i) = field(X(i,1),X(i,2));
end

VAR = variogram([X(:,1), X(:,2)], Z, 'plot', false, 'type', 'gamma', 'nrbins', ceil(length(Z)/2));
c0 = max(VAR.val);
a0 = max(VAR.distance)*2/3;

ix = linspace(0, max(VAR.distance));

kernels = {'gaussian'};

figure()
plots = zeros(1, max(size(kernels)) + 1);

for i = 1:max(size(kernels))
    [range, sill, nugget, S] = variogramfit(VAR.distance, VAR.val, a0, c0, [], 'model', char(kernels(i)), 'plotit', false);
    b = [range, sill, nugget];
    kern = S.func(b, ix);
    plots(i) = plot(ix, kern, 'k-', 'LineWidth', 2);
    hold on
end

plot(VAR.distance, VAR.val, 'ro', 'Linewidth', 1.1);

title('Variogram Fit Semi-Variogram')
xlabel('Distance')
ylabel('\gamma', 'FontSize', 17)

axis([0, max(VAR.distance), 0, max(VAR.val)])


hold on
plot([0 range], [sill sill], '--k')
hold on
plot([range range], [0 sill], '--r')

legend('Gaussian Model', 'Semi-Variogram', 'range', 'sill', 'northwest')

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

warning('off','all');

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
        
        idx = knnsearch(X, [i, j], 'k', 1);
        neighbor_conf(idx, 1) = sum(d);
        neighbor_conf(idx, 2) = neighbor_conf(idx, 2) + 1;
    end
end

for i = 1:length(neighbor_conf)
    if (neighbor_conf(i, 2) == 0)
        %neighbor_conf(i, 1) = neighbor_conf(i, 1);
        neighbor_conf(i, 1) = 1;
    else
        %neighbor_conf(i, 1) = neighbor_conf(i, 1) / neighbor_conf(i, 2); % calculate average uncertainty of each neighbor
        neighbor_conf(i, 1) = 1 / neighbor_conf(i, 2);
    end
end


figure()

pcolor(pred_field)
title('Kriging Predicted Field')
shading interp; % gets rid of the grid lines on the surf()
xlabel('x_1');
ylabel('x_2');
zlabel('u');
%%
hold on

[vx,vy] = voronoi(X(:,1),X(:,2));

% Assign labels to the points.
nump = length(X);
plabels = arrayfun(@(n) {sprintf('N_{%d}\nG_{%d}=%f', n,n,neighbor_conf(n,1))}, (1:nump)');
Hpl = text(X(:,1), X(:,2), plabels, 'FontWeight', ...
      'bold', 'HorizontalAlignment','center');

title('Natural Neighbor Selection of the predicted field');
plot(vx, vy, 'k-', 'LineWidth', 1.1)

%%
nn = 5;

edges = [];
weights = [];
for curr_n = 1:length(X)

    ix = knnsearch(X,X(curr_n,:),'k', nn + 1);
    for j = 2:length(ix)

        if (length(edges) > 0)
            if ( (sum(ismember(edges, [curr_n, ix(j)], 'rows')) == 0) && (sum(ismember(edges, [ix(j), curr_n], 'rows')) == 0) )
                edges = [edges; [curr_n, ix(j)]];
                weights = [weights ; 1/(abs(neighbor_conf(curr_n,1)) + abs(neighbor_conf(ix(j),1))) ];
            end
        else
            edges = [edges; [curr_n, ix(j)]];
            weights = [weights ; 1/(abs(neighbor_conf(curr_n,1)) + abs(neighbor_conf(ix(j),1))) ];
        end
    end
end

G = graph(edges(:,1)', edges(:,2)', weights');
figure()

Gp = plot(G,'ko', 'LineWidth', 2);%, 'EdgeLabel', G.Edges.Weight)
Gp.Marker = 'o';
Gp.NodeColor = 'r';
Gp.MarkerSize = 9;
% turn off graph axes
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Directed Graph Representation of Voronoi Tesselated Field')

%%

error = 0;

for i = 1:length(field)
    for j = 1:length(field)
        error = error + (field(i,j) - pred_field(i,j))^2;
    end
end
error = error / length(field)^2;
fprintf('For %d Samples, %d iterations:\n', n, runs)
fprintf('Error_RMS = %f\n', error)

errors(runs) = error;
uncerts(runs) = mean(neighbor_conf(:,1));
sample_track(runs) = n;
%%

current_neighbor = knnsearch(X,X(length(X),:),'k',1);
worst_neighbor = find(neighbor_conf == min(neighbor_conf(:,1)));
uav_path = shortestpath(G,current_neighbor(1),worst_neighbor(1));

new_samples = zeros(length(uav_path), 2);

% for i = 1:length(uav_path)
%     new_samples(i,:) = X(uav_path(i),:) + randi(20, [1 2]);
%     if (new_samples(i,1) > length(field))
%         new_samples(i,1) = length(field);
%     end
%     if (new_samples(i,2) > length(field))
%         new_samples(i,2) = length(field);
%     end
% end

n = n + 1;
X = [randi(length(field), n+length(uav_path), 1), randi(length(field), n+length(uav_path), 1)];

end
%%
figure()
plot(1:length(errors), errors, 'r', 'LineWidth', 2); 
axis([1 length(errors) 0 max(errors)]); 
title('Error Over Iterations',  'Interpreter', 'latex', 'FontSize', 20); 
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Error', 'Interpreter', 'latex', 'FontSize', 20); 
hold on; 
p = polyfit(1:length(errors), errors', 1); 
plot(1:length(errors),polyval(p, 1:length(errors))', ':k', 'LineWidth', 2);
grid on;
export_img_latex(gcf, '../figures/error_plot')

save('good_err.mat', 'errors')
