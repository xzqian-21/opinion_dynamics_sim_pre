%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig1b.m
%
% This script generates Fig.1(b) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
% rng(43);

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

alpha_arr  = [1/6, 1/2, 5/6];

output_dir = './fig';
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

    % --- Opinion pattern ---
    PatternPlot(B, n);  % use full 50x50 lattice
    opinion_filename = fullfile(output_dir, sprintf('fig1_c%d.svg', idx));
    saveas(gcf, opinion_filename);

    % --- Conflict pattern ---
    ConflictPlot(B, n);
    conflict_filename = fullfile(output_dir, sprintf('fig1_d%d.svg', idx));
    saveas(gcf, conflict_filename);
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

