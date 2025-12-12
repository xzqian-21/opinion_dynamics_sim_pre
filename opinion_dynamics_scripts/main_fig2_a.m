%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig2a.m
%
% This script generates Fig.2(a) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  Simulation: compute disagreement & polarization metrics over alpha
% -------------------------------------------------------------------------
clear; clc; close all;

% Network size
N = 2500;           % total number of nodes (50x50 lattice)
n = 50;             % lattice side length 

data_file = './network_matrix/R50.mat';
S = load(data_file);
A = S.A;

deg = 4 * ones(N,1);                     
W   = A / 4;                              
L   = sparse(spdiags(deg,0,N,N) - A);     

alpha_grid = (0.5:0.1:2.5) / 3;           
n_alpha    = numel(alpha_grid);

% Storage for Monte Carlo repetitions
n_rep   = 100;
rec_2d  = zeros(n_rep, n_alpha);          
rec_2d1 = zeros(n_rep, n_alpha);          

rng(43);                                

tic;
parfor j = 1:n_rep
    [a, a1] = SimuScr(A, W, L, alpha_grid);
    rec_2d(j,:)  = a;
    rec_2d1(j,:) = a1;
end
toc;

% Save all relevant variables for plotting

save('./result/2d.mat', 'rec_2d', 'rec_2d1', 'alpha_grid');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

load('./result/2d.mat');   % rec_2d, rec_2d1, alpha_grid
load('color.mat');               % expects: color1, color3, color12, color32

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
tick_size  = 18;
label_size = 18;

hold on;

% Vertical reference lines at alpha = 1/6 and 5/6
l1 = plot([1/6, 1/6], [0, 1], 'LineWidth', 4, 'LineStyle', ':');
l2 = plot([5/6, 5/6], [0, 1], 'LineWidth', 4, 'LineStyle', ':');
l1.Color = '#8a2be2';
l2.Color = '#8a2be2';

% We sweep alpha from high to low in the original code; to keep the same
% direction, reverse along the alpha dimension
rec_2d_plot  = rec_2d(:, end:-1:1);
rec_2d1_plot = rec_2d1(:, end:-1:1);

alpha_plot = alpha_grid(end:-1:1);

% Disagreement metric: mean ± std
mu_dis  = mean(rec_2d_plot, 1);
std_dis = std(rec_2d_plot, 0, 1);

e = errorbar(alpha_plot, mu_dis, std_dis, std_dis, ...
    'DisplayName', 'Disagree. Metric');
e.Color           = color3;
e.LineStyle       = '-';
e.LineWidth       = 4/3;
e.Marker          = 'square';
e.MarkerSize      = 10;
e.MarkerEdgeColor = color32;
e.MarkerFaceColor = color3;

% Polarization metric: mean ± std
mu_pol  = mean(rec_2d1_plot, 1);
std_pol = std(rec_2d1_plot, 0, 1);

e1 = errorbar(alpha_plot, mu_pol, std_pol, std_pol, ...
    'DisplayName', 'Polarization Metric');
e1.Color           = color1;
e1.LineStyle       = '-';
e1.LineWidth       = 4/3;
e1.Marker          = 'square';
e1.MarkerSize      = 10;
e1.MarkerEdgeColor = color12;
e1.MarkerFaceColor = color1;

% Axes limits and ticks
ylim([0, 1]);
yticks(0:0.25:1);
yticklabels({'0','','0.5','','1'});

xlim([0.1, 0.9]);
xticks([0.2, 0.4, 0.6, 0.8, 0.9]);
xticklabels({'0.2','0.4','0.6','0.8',''});

grid on;
box off;
ax = gca;

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(gca, 'Color', 'none');

legend([e, e1], 'FontSize', 18, 'FontName', 'Arial', 'Location', 'Northwest');

xlabel('Strength {\alpha}', 'FontName', 'Arial', 'FontSize', label_size);
% ylabel('Index Value','FontName','Arial','FontSize',label_size);

if ~exist('./fig','dir')
    mkdir('./fig');
end
saveas(gcf, './fig/fig2_a.svg');

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
    y_new  = Y(end,:).';              % last time point

    res = max(abs(y_new - y));        % max-norm distance

    y = y_new;
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
