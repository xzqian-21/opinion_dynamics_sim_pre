%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig4_f.m
%
% This script generates Fig.4(f) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  Simulation: collect spectra before and after selective node removal
% -------------------------------------------------------------------------
clear; clc; close all;

N         = 2500;
m         = 5;
T         = 0.5;
gamma     = 3.1;
distr     = 0;
plot_flag = 0;

n_total   = 100;
n_remove  = 100;

rec_sp_old = zeros(n_total, N);  
rec_sp_new = -ones(n_total, N);  

tic;
parfor j = 1:n_total
    rng(j + 100);

    
    [A, coords, comm, d] = nPSO_model(N, m, T, gamma, distr, plot_flag); 
    A   = sparse(A);
    deg = A * ones(N,1);

    W     = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L_bar = speye(N) - W;

    sp_old           = eig(full(L_bar));
    rec_sp_old(j,:)  = sp_old;

    c          = 2.5;          
    d_nbr      = 0.5;         
    K          = d_nbr * W;

    mean_kappa = 3.0;
    sigma      = 0.2;
    var_kappa  = randn(N,1);
    kappa      = mean_kappa*(1+ sigma * var_kappa);

    E0 = randn(N,1);
    y0 = [E0 / 5];

    B = SNevo(kappa, c, K, y0);
    loc_conf = (A .* (B' - B).^2) * ones(N,1);
    [~, mask] = maxk(loc_conf, n_remove);  

    A_new = full(A);
    A_new(mask, :) = [];
    A_new(:, mask) = [];

    d_vec = A_new * ones(N - n_remove, 1);
    D_new = diag(d_vec);

   
    zero_idx = find(d_vec == 0);
    if ~isempty(zero_idx)
        D_new(zero_idx, :) = [];
        D_new(:, zero_idx) = [];
        A_new(zero_idx, :) = [];
        A_new(:,zero_idx) = [];
    end

    N_red = size(D_new, 1);
    if N_red == 0
        continue;
    end

    W_new     = D_new^(-1) * A_new;
    L_bar_new = eye(N_red) - W_new;

    sp_new = eig(L_bar_new);
    sp_new_padded   = [sp_new; -ones(N - N_red, 1)];
    rec_sp_new(j,:) = sp_new_padded;
end
toc;

rec_sp_new = rec_sp_new(:);
rec_sp_new = real(rec_sp_new);
rec_sp_new = sort(rec_sp_new);

rec_sp_old = real(rec_sp_old(:));

save('./result/sp_select_100.mat', ...
     'rec_sp_new', 'rec_sp_old', 'n_total', 'N', 'n_remove');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot: before vs after eigenvalue distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
load('./result/sp_select_100.mat');   
load('color.mat');

% -------------------------------------------------------------------------
%  Build histograms and smoothed densities
% -------------------------------------------------------------------------
bin_edges   = 0:0.01:2;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

valid_new = rec_sp_new(rec_sp_new > -0.5);

ratio = numel(rec_sp_old) / numel(valid_new);

[counts_new, ~] = histcounts(valid_new,  bin_edges);
[counts_old, ~] = histcounts(rec_sp_old, bin_edges);

height1 = movmean(counts_new, 4) * ratio;  % after removal
height2 = movmean(counts_old, 4);          % before removal

idx_patch = find(bin_centers < 1.5);
x_patch   = bin_centers(idx_patch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Combined figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
ax = axes('Parent', gcf);
ax.Position = [0.2, 0.25, 0.6, 0.6];

tick_size  = 18;
label_size = 18;
hold(ax, 'on');

% Filled areas (before & after)
patch([x_patch, fliplr(x_patch)], ...
      [height2(idx_patch), zeros(1, numel(idx_patch))], ...
      color3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

patch([x_patch, fliplr(x_patch)], ...
      [height1(idx_patch), zeros(1, numel(idx_patch))], ...
      color1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Smoothed curves
l2 = plot(ax, bin_centers, height2, ...
    'LineWidth', 4, 'Color', color32, 'DisplayName', 'Before');
l1 = plot(ax, bin_centers, height1, ...
    'LineWidth', 4, 'Color', color12, 'DisplayName', 'After');

grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, 1]);
ylim(ax, [0, 37.5 * n_total]);

xticks(ax, 0:0.25:1);
yticks(ax, 0:12.5*n_total:37.5*n_total);
yticklabels(ax, {'0','0.5','1','1.5'});

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);

xlabel('Eigenvalue {\lambda}', 'FontName', 'Arial', 'FontSize', label_size);
ylabel('Empirical PDF',       'FontName', 'Arial', 'FontSize', label_size);

legend([l2, l1], 'FontName', 'Arial', 'FontSize', tick_size, ...
       'Location', 'Northwest');

set(ax, 'Color', 'none');

saveas(gcf, './fig/fig4_f.svg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dydt = SN(t, y, kappa, c, K)


B = y(:);
N = numel(kappa);
if numel(B) ~= N
    error('Length of y (%d) does not match length of kappa (%d).', numel(B), N);
end

F_ext = 0;
E     = c .* B + K * B + F_ext;
dydt  = -B + 2 ./ (1 + exp(-kappa .* E)) - 1;
end


function B = SNevo(kappa, c, K, y0, tspan, tol, max_iter)


N = numel(kappa);

if nargin < 5 || isempty(tspan)
    tspan = 0:5:250;
end
if nargin < 6 || isempty(tol)
    tol = 1e-3;
end
if nargin < 7 || isempty(max_iter)
    max_iter = 512;
end

opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);

y = y0(:);
if numel(y) ~= N
    error('Length of y0 (%d) does not match length of kappa (%d).', numel(y), N);
end

for count = 1:max_iter
    [~, Y] = ode45(@(t, y) SN(t, y, kappa, c, K), tspan, y, opts);
    y_new  = Y(end,:).';

    res = max(abs(y_new - y));
    y   = y_new;
    if res < tol
        break;
    end
end

B = y(1:N);
end
