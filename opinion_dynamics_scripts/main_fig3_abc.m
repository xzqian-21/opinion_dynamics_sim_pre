%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig3_abc.m
%
% This script generates Fig.3(a)-(c) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------
clear; clc; close all;

N        = 2500;
alpha_grid = (0.5:0.1:2.5) / 3;   % 21 points from ~1/6 to ~5/6
n_alpha  = numel(alpha_grid);
n_rep    = 100;


%% ------------------------------------------------------------------------
%  Experiment 1: Watts–Strogatz small-world network
% -------------------------------------------------------------------------
rec_sw  = zeros(n_rep, n_alpha);   % conflict metric
rec_sw1 = zeros(n_rep, n_alpha);   % polarization metric

tic;
parfor j = 1:n_rep
    rng(43+j)
    A = sw_net('N', N, 'k', 4, 'p', 0.1);
    A = sparse(A);
    deg = A * ones(N,1);

    W = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L = sparse(spdiags(deg, 0, N, N) - A);

    [a, a1] = SimuScr(A, W, L, alpha_grid);
    rec_sw(j,:)  = a;
    rec_sw1(j,:) = a1;
end
toc;
save('./result/sw.mat', 'rec_sw', 'rec_sw1', 'alpha_grid');


%% ------------------------------------------------------------------------
%  Experiment 2: Barabási–Albert scale-free network
% -------------------------------------------------------------------------
clearvars -except N alpha_grid n_alpha n_rep;


rec_ba  = zeros(n_rep, n_alpha);
rec_ba1 = zeros(n_rep, n_alpha);

tic;
parfor j = 1:n_rep
    rng(43+j)
    A = ba_net('N', N, 'm0', 2, 'm', 2);
    A = sparse(A);
    deg = A * ones(N,1);

    W = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L = sparse(spdiags(deg, 0, N, N) - A);

    [a, a1] = SimuScr(A, W, L, alpha_grid);
    rec_ba(j,:)  = a;
    rec_ba1(j,:) = a1;
end
toc;
save('./result/ba.mat', 'rec_ba', 'rec_ba1', 'alpha_grid');


%% ------------------------------------------------------------------------
%  Experiment 3: PSO / hyperbolic network
% -------------------------------------------------------------------------
clearvars -except N alpha_grid n_alpha n_rep;
rng(43);

rec_ps  = zeros(n_rep, n_alpha);
rec_ps1 = zeros(n_rep, n_alpha);

m         = 2;
T         = 0.2;
gamma     = 2.1;
distr     = 0;
plot_flag = 0;

tic;
parfor j = 1:n_rep

    rng(43+j)
    [A, coords, comm, d] = nPSO_model(N, m, T, gamma, distr, plot_flag); 
    A = sparse(A);
    deg = A * ones(N,1);

    W = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L = sparse(spdiags(deg, 0, N, N) - A);

    [a, a1] = SimuScr(A, W, L, alpha_grid);
    rec_ps(j,:)  = a;
    rec_ps1(j,:) = a1;
end
toc;
save('./result/ps.mat', 'rec_ps', 'rec_ps1', 'alpha_grid');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot: WS network (Fig.4 – WS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
load('./result/sw.mat');   
load('color.mat');            

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
tick_size  = 18;
label_size = 18;

ax = axes('Parent', gcf);
ax.Position = [0.2, 0.25, 0.6, 0.6];

hold(ax, 'on');

alpha_plot  = alpha_grid(end:-1:1);
rec_sw_plot  = rec_sw(:, end:-1:1);
rec_sw1_plot = rec_sw1(:, end:-1:1);

mu_conf  = mean(rec_sw_plot, 1);
std_conf = std(rec_sw_plot, 0, 1);

e = errorbar(alpha_plot, mu_conf, std_conf, std_conf, ...
    'DisplayName', 'Disag. Metric');
e.Color           = color3;
e.LineStyle       = '-';
e.LineWidth       = 4/3;
e.Marker          = 'square';
e.MarkerSize      = 10;
e.MarkerEdgeColor = color32;
e.MarkerFaceColor = color3;

mu_pol  = mean(rec_sw1_plot, 1);
std_pol = std(rec_sw1_plot, 0, 1);
alpha_pol = alpha_plot + 0.02/3;   % corresponds to (.52:0.1:2.52)/3

e1 = errorbar(alpha_pol, mu_pol, std_pol, std_pol, ...
    'DisplayName', 'Polar. Metric');
e1.Color           = color1;
e1.LineStyle       = '-';
e1.LineWidth       = 4/3;
e1.Marker          = 'square';
e1.MarkerSize      = 10;
e1.MarkerEdgeColor = color12;
e1.MarkerFaceColor = color1;

ylim([0, 1.0]);
yticks(0:0.25:1);
yticklabels({'0','','0.5','','1'});

xlim([0.1, 0.9]);
xticks([0.2:0.2:0.8, 0.9]);
xticklabels({'0.2','0.4','0.6','0.8',''});

grid on; box off;

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(gca, 'Color', 'none');

xlabel('Strength {\alpha}', 'FontName', 'Arial', 'FontSize', label_size);

saveas(gcf, './fig/fig3_ws.svg');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot: BA network (Fig.4 – BA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
load('./result/ba.mat');   % rec_ba, rec_ba1, alpha_grid
load('color.mat');

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
tick_size  = 18;
label_size = 18;

ax = axes('Parent', gcf);
ax.Position = [0.2, 0.25, 0.6, 0.6];

rec_ba_plot  = rec_ba(:, end:-1:1);
rec_ba1_plot = rec_ba1(:, end:-1:1);
alpha_plot   = alpha_grid(end:-1:1);

% Conflict index
mu_conf  = mean(rec_ba_plot, 1);
std_conf = std(rec_ba_plot, 0, 1);

e = errorbar(alpha_plot, mu_conf, std_conf, std_conf, ...
    'DisplayName', 'Disag. Metric');
e.Color           = color3;
e.LineStyle       = '-';
e.LineWidth       = 4/3;
e.Marker          = 'square';
e.MarkerSize      = 10;
e.MarkerEdgeColor = color32;
e.MarkerFaceColor = color3;

hold(ax, 'on');

% Polarization index
mu_pol  = mean(rec_ba1_plot, 1);
std_pol = std(rec_ba1_plot, 0, 1);
alpha_pol = alpha_plot + 0.02/3;

e1 = errorbar(alpha_pol, mu_pol, std_pol, std_pol, ...
    'DisplayName', 'Polar. Metric');
e1.Color           = color1;
e1.LineStyle       = '-';
e1.LineWidth       = 4/3;
e1.Marker          = 'square';
e1.MarkerSize      = 10;
e1.MarkerEdgeColor = color12;
e1.MarkerFaceColor = color1;

ylim([0, 1.0]);
yticks(0:0.25:1);
yticklabels({'0','','0.5','','1'});

xlim([0.1, 0.9]);
xticks([0.2:0.2:0.8, 0.9]);
xticklabels({'0.2','0.4','0.6','0.8',''});

grid on; box off;

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(gca, 'Color', 'none');

xlabel('Strength {\alpha}', 'FontName', 'Arial', 'FontSize', label_size);

saveas(gcf, './fig/fig3_ba.svg');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot: PS network (Fig.4 – PS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
load('./result/ps.mat');   % rec_ps, rec_ps1, alpha_grid
load('color.mat');

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
tick_size  = 18;
label_size = 18;

ax = axes('Parent', gcf);
ax.Position = [0.2, 0.25, 0.6, 0.6];

rec_ps_plot  = rec_ps(:, end:-1:1);
rec_ps1_plot = rec_ps1(:, end:-1:1);
alpha_plot   = alpha_grid(end:-1:1);

% Conflict index
mu_conf  = mean(rec_ps_plot, 1);
std_conf = std(rec_ps_plot, 0, 1);

e = errorbar(alpha_plot, mu_conf, std_conf, std_conf, ...
    'DisplayName', 'Disag. Metric');
e.Color           = color3;
e.LineStyle       = '-';
e.LineWidth       = 4/3;
e.Marker          = 'square';
e.MarkerSize      = 10;
e.MarkerEdgeColor = color32;
e.MarkerFaceColor = color3;

hold(ax, 'on');

% Polarization index
mu_pol  = mean(rec_ps1_plot, 1);
std_pol = std(rec_ps1_plot, 0, 1);
alpha_pol = alpha_plot + 0.02/3;

e1 = errorbar(alpha_pol, mu_pol, std_pol, std_pol, ...
    'DisplayName', 'Polar. Metric');
e1.Color           = color1;
e1.LineStyle       = '-';
e1.LineWidth       = 4/3;
e1.Marker          = 'square';
e1.MarkerSize      = 10;
e1.MarkerEdgeColor = color12;
e1.MarkerFaceColor = color1;

ylim([0, 1.0]);
yticks(0:0.25:1);
yticklabels({'0','','0.5','','1'});

xlim([0.1, 0.9]);
xticks([0.2:0.2:0.8, 0.9]);
xticklabels({'0.2','0.4','0.6','0.8',''});

grid on; box off;

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(gca, 'Color', 'none');

xlabel('Strength {\alpha}', 'FontName', 'Arial', 'FontSize', label_size);
saveas(gcf, './fig/fig3_ps.svg');


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