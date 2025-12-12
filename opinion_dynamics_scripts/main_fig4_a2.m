%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig4_a2.m
%
% This script generates Fig.4(a2) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------

clear; clc; close all;


% Alpha grid (21 points, from ~1/6 to ~5/6)
alpha_grid = (0.5:0.1:2.5) / 3;
n_alpha    = numel(alpha_grid);


n_rep  = 100;
rec_tm  = zeros(n_rep, n_alpha);   
rec_tm1 = zeros(n_rep, n_alpha);  

tic;
parfor j = 1:n_rep

    rng(j + 100);

    N         = 2500;
    m         = 5;
    T         = 0.5;
    gamma     = 3.1;
    distr     = 0;
    plot_flag = 0;

    [A, coords, comm, d] = nPSO_model(N, m, T, gamma, distr, plot_flag); 
    A   = sparse(A);
    deg = A * ones(N,1); 
 
    n_remove = 100;
    A_new    = A(n_remove+1:N, n_remove+1:N);
    d_vec    = A_new * ones(N - n_remove, 1);
    zero_idx = find(d_vec == 0);

    if ~isempty(zero_idx)
        keep_idx = setdiff(1:(N - n_remove), zero_idx);
        A_new    = A_new(keep_idx, keep_idx);
    end

    N_red = size(A_new, 1);
    d_vec = A_new * ones(N_red, 1);
    D_new = diag(d_vec);

    W_new = D_new^(-1) * A_new;
    L_new = D_new - A_new;

    A_red = sparse(A_new);
    W_red = sparse(W_new);
    L_red = sparse(L_new);

    [a, a1] = SimuScr(A_red, W_red, L_red, alpha_grid);

    rec_tm(j,:)  = a;
    rec_tm1(j,:) = a1;
end
toc;

% Reverse alpha order to match原来从大到小的绘图习惯
rec_tm  = rec_tm(:, end:-1:1);
rec_tm1 = rec_tm1(:, end:-1:1);

save('./result/CI_str_100.mat', 'rec_tm', 'rec_tm1', 'alpha_grid', 'n_rep');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot C and P vs alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
load('./result/CI_str_100.mat');  % rec_tm, rec_tm1, alpha_grid, n_rep
load('color.mat');

figure('OuterPosition',[800, 200, 320, 350], 'visible', 'on');
tick_size  = 18;
label_size = 18;

ax = axes('Parent', gcf);
ax.Position = [0.2, 0.25, 0.6, 0.6];

rec_tm_sub  = rec_tm(:, 1:2:21);
rec_tm1_sub = rec_tm1(:, 1:2:21);

alpha_full = alpha_grid(end:-1:1);      
alpha_sub  = alpha_full(1:2:21);       
alpha_ci   = alpha_sub;                
alpha_pi   = alpha_sub + 0.02/3;       

% CI: mean ± std
mu_ci  = mean(rec_tm_sub, 1);
std_ci = std(rec_tm_sub, 0, 1);

hold(ax, 'on');
e = errorbar(ax, alpha_ci, mu_ci, std_ci, std_ci, ...
    'DisplayName', 'Disag. Metric');
e.Color           = color3;
e.LineStyle       = '-';
e.LineWidth       = 4/3;
e.Marker          = 'square';
e.MarkerSize      = 10;
e.MarkerEdgeColor = color32;
e.MarkerFaceColor = color3;

% PI: mean ± std
mu_pi  = mean(rec_tm1_sub, 1);
std_pi = std(rec_tm1_sub, 0, 1);

e1 = errorbar(ax, alpha_pi, mu_pi, std_pi, std_pi, ...
    'DisplayName', 'Polar. Metric');
e1.Color           = color1;
e1.LineStyle       = '-';
e1.LineWidth       = 4/3;
e1.Marker          = 'square';
e1.MarkerSize      = 10;
e1.MarkerEdgeColor = color12;
e1.MarkerFaceColor = color1;

% Axis settings
ylim(ax, [0, 1]);
yticks(ax, 0:0.25:1);
yticklabels(ax, {'0','','0.5','','1'});

xlim(ax, [0.1, 0.9]);
xticks(ax, [0.1, 0.2:0.2:0.8, 0.9]);
xticklabels(ax, {'', '0.2','0.4','0.6','0.8', ''});

grid(ax, 'on'); box(ax, 'off');

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(gca, 'Color', 'none');

legend('FontSize', 18, 'FontName', 'Arial', 'Location', 'Northwest');
xlabel('Strength {\alpha}', 'FontName', 'Arial', 'FontSize', label_size);

saveas(gcf, './fig/fig4_a2.svg');

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


function [record_arr, record_arr1] = SimuScr(A, W, L, alpha_grid)


N        = size(A,1);
n_alpha  = numel(alpha_grid);

record_arr  = zeros(1, n_alpha);
record_arr1 = zeros(1, n_alpha);

B0 = zeros(N,1);
E0 = randn(N,1);
y0 = E0 / 5;

mean_kappa = 3;
sigma      = 0.2;
var_kappa  = randn(N,1);
kappa      = mean_kappa*(1 + sigma * var_kappa);

for i = 1:n_alpha
    alpha = alpha_grid(i);

    c = alpha;
    d = 1 - alpha;
    K = d * W;

    B = SNevo(kappa, c, K, y0);

    record_arr(i)  = B' * L * B / sum(A(:));
    record_arr1(i) = var(B(:));
end
end
