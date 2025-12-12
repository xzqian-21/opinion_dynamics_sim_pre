%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig1b.m
%
% This script generates Fig.2(b)and(c) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
rng(43);

%% ------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------
N = 2500;              % total number of nodes (50x50 lattice)
n = 50;                % lattice side length

data_file = './network_matrix/R50.mat';
S = load(data_file);   
A = S.A;
W = A / 4;             

mean_kappa = 3;
sigma      = 0.2;
var_kappa  = randn(N,1);                      
kappa      = mean_kappa * (1 + sigma*var_kappa);

alpha_arr  = [1/6, 5/6];

for idx = 1:numel(alpha_arr)
    alpha = alpha_arr(idx);

    % Local field parameters
    c = alpha;
    d = (1 - alpha);
    K = d * W;  % effective coupling matrix

    % Initial condition
    B0 = randn(N,1);
    y0 = B0 / 5;

    % Time evolution to steady state
    B = SNevo(kappa, c, K, y0);
    DistributionPlot(B, W);

    if idx == 1
        saveas(gcf, './fig/fig2_b.svg');
    elseif idx == 2
        saveas(gcf, './fig/fig2_c.svg');
    end
end

disp('Simulation and figure generation finished.');

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
E = c .* B + K * B + F_ext;
dydt = -B + 2 ./ (1 + exp(-kappa .* E)) - 1;
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

opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

y = y0(:);
if numel(y) ~= N
    error('Length of y0 (%d) does not match length of kappa (%d).', numel(y), N);
end

for count = 1:max_iter
    [~, Y] = ode45(@(t, y) SN(t, y, kappa, c, K), tspan, y, opts);
    y_new  = Y(end,:).';  

    res = max(abs(y_new - y)); 

    y = y_new;
    if res < tol
        break;
    end
end

B = y(1:N);
end

function PatternPlot(B, n, ns)

if nargin < 3 || isempty(ns)
    ns = n;
end
ns = min(ns, n);

B = B(:);
if numel(B) ~= n^2
    error('Length of B (%d) does not match n^2 (%d).', numel(B), n^2);
end

B_grid = reshape(B, n, n);
B_sub  = B_grid(1:ns, 1:ns);

figure('OuterPosition',[800, 200, 480, 480], 'visible', 'on');
imagesc(1:ns, 1:ns, B_sub.');
axis square;
hold on;

C     = (-1:0.01:1).';
b_end = [57,  91, 139] / 255;   
w_end = [245, 245, 245] / 255;   
y_end = [245, 224, 166] / 255;   

c_l = (C > 0) .* C      * b_end + ...
      (C < 0) .* (-C)   * y_end + ...
      (1 - abs(C))      * w_end;

colormap(c_l);
clim([-1, 1]);

set(gca, 'FontSize', 24, 'color', 'none', 'linewidth', 1);
xticks([]); yticks([]);
xlabel('', 'FontSize', 24, 'FontName', 'Arial');
ylabel('', 'FontSize', 24, 'FontName', 'Arial');
axis square;
end

function ConflictPlot(B, n)

B = B(:);
if numel(B) ~= n^2
    error('Length of B (%d) does not match n^2 (%d).', numel(B), n^2);
end

B_grid = reshape(B, n, n);

B_l = B_grid - circshift(B_grid, [0, -1]);   
B_r = B_grid - circshift(B_grid, [0,  1]);   
B_u = B_grid - circshift(B_grid, [-1, 0]);  
B_d = B_grid - circshift(B_grid, [ 1, 0]);   


B_c = abs(B_l.^2) + abs(B_r.^2) + abs(B_u.^2) + abs(B_d.^2);

figure('OuterPosition',[800, 200, 460, 460], 'visible', 'on');
imagesc(1:n, 1:n, B_c.');
axis square;
hold on;

set(gca, 'FontSize', 24, 'color', 'none', 'linewidth', 1);
xticks([]); yticks([]);

interval = (0:0.01:1).';
s = [255, 255, 255] / 255;   
e = [255,   0,   0] / 255;  
RGB = interval * e + (1 - interval) * s;
colormap(RGB);
clim([0, 4]);

xlabel('', 'FontSize', 24, 'FontName', 'Arial');
ylabel('', 'FontSize', 24, 'FontName', 'Arial');
end

function DistributionPlot(B, W)
% DistributionPlot  Opinion vs neighborhood mean with marginal densities.
%
%   DistributionPlot(B, W)
%
%   Inputs:
%       B  - N x 1 vector of final opinions
%       W  - N x N weight matrix; W*B is interpreted as local mean opinion

B = B(:);
N = numel(B);

if size(W,1) ~= N || size(W,2) ~= N
    error('W must be an N x N matrix with N = length(B).');
end

if ~isfile('color.mat')
    error('File "color.mat" not found. It must define color1, color3, color12, color32.');
end
load('color.mat');  % expects color1, color3, color12, color32

x  = B;
yN = W * B;   % neighborhood mean opinion

figure('OuterPosition',[800, 200, 380, 460], 'visible', 'on'); % original 420x500

%% Main scatter plot
ax1 = axes('Parent', gcf);
hold(ax1, 'on');

scatter(ax1, x, yN, 20, 'filled', ...
    'MarkerEdgeColor',  color1, ...
    'MarkerFaceColor',  color1, ...
    'MarkerEdgeAlpha',  0.3, ...
    'MarkerFaceAlpha',  0.3);

ax1.Position = [0.2, 0.2, 0.6, 0.6];
ax1.XLim     = [-1, 1];
ax1.YLim     = [-1, 1];
ax1.XTick    = -1:0.5:1;
ax1.YTick    = -1:0.5:1;
ax1.XTickLabel = {'-1','','0','','1'};
ax1.YTickLabel = {'-1','','0','','1'};

tick_size = 24;
font_size = 24;

set(ax1, 'box', 'on', 'XGrid', 'on', 'YGrid', 'on');
set(ax1.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax1.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);

xlabel(ax1, '$\boldmath{x}_i$',              'Interpreter', 'latex', 'FontSize', font_size);
ylabel(ax1, '$\boldmath{\bar{x}}_{\mathcal{N}_i}$', 'Interpreter', 'latex', 'FontSize', font_size);

%% Top marginal density of B
ax2 = axes('Parent', gcf);
hold(ax2, 'on');

[f, xi] = ksdensity(x, 'Bandwidth', 0.05);
fill(ax2, [xi, xi(1)], [f, 0], color1, ...
    'FaceAlpha', 0.3, 'EdgeColor', color12, 'LineWidth', 1.2);

ax2.Position   = [0.2, 0.82, 0.6, 0.1];
ax2.YColor     = 'none';
ax2.XTickLabel = '';
ax2.TickDir    = 'out';
ax2.XLim       = ax1.XLim;
ax2.XTick      = ax1.XTick;

%% Right marginal density of W*B
ax3 = axes('Parent', gcf);
hold(ax3, 'on');

[f2, yi] = ksdensity(yN, 'Bandwidth', 0.05);
fill(ax3, [f2, 0], [yi, yi(1)], color3, ...
    'FaceAlpha', 0.3, 'EdgeColor', color32, 'LineWidth', 1.2);

ax3.Position   = [0.82, 0.2, 0.1, 0.6];
ax3.XColor     = 'none';
ax3.YTickLabel = '';
ax3.TickDir    = 'out';
ax3.YLim       = ax1.YLim;
ax3.YTick      = ax1.YTick;

set(ax1, 'color', 'none', 'linewidth', 4/3);
set(ax2, 'color', 'none');
set(ax3, 'color', 'none');
end

