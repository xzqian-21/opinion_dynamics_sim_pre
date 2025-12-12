%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig4_c.m
%
% This script generates Fig.4(c) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  Simulation: CI under node removal (By Order vs By Disag.)
% -------------------------------------------------------------------------
clear; clc; close all;


n_rep = 100;          
n_steps = 11;         % 0, 10, 20, ..., 100 nodes removed

recA_conf = zeros(n_rep, n_steps); 
recA_deg  = zeros(n_rep, n_steps);  
recB_conf = zeros(n_rep, n_steps); 
recB_deg  = zeros(n_rep, n_steps);  


N         = 2500;
m         = 5;
T         = 0.5;
gamma     = 3.1;
distr     = 0;
plot_flag = 0;   
alpha      = 5/6;      
c          = alpha;
d_nbr      = 1 - alpha;
mean_kappa = 3;
sigma      = 0.2;

tic;
parfor i = 1:n_rep
    rng(i + 100);

  
    [A, coords, comm, d_param] = nPSO_model(N, m, T, gamma, distr, plot_flag); 
    A   = sparse(A);
    deg = A * ones(N,1);
    D   = diag(deg);

    W = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L = sparse(D - A);

    
    var_kappa = randn(N,1);
    kappa     = mean_kappa * (1 + sigma * var_kappa);

    K = d_nbr * W;

    B0 = randn(N,1);
    y0 = B0 / 5;

    B_full = SNevo(kappa, c, K, y0);

    ci0   = B_full' * L * B_full / sum(A(:));     
    deg0  = sum(A(:)) / N;                        

    recA_conf_i = zeros(1, n_steps);
    recA_deg_i  = zeros(1, n_steps);
    recB_conf_i = zeros(1, n_steps);
    recB_deg_i  = zeros(1, n_steps);

    recA_conf_i(1) = ci0;
    recA_deg_i(1)  = deg0;
    recB_conf_i(1) = ci0;
    recB_deg_i(1)  = deg0;

    for j = 1:10
        n_rm = 10 * j;

        [A_new, D_new, B_new, kappa_new] = params_filter_A(n_rm, A, D, B_full, kappa);

        N_red = size(D_new,1);
        if N_red == 0
            recA_conf_i(j+1) = NaN;
            recA_deg_i(j+1)  = NaN;
            continue;
        end

        W_new = D_new^(-1) * A_new;
        L_new = D_new - A_new;

        W_new = sparse(W_new);
        L_new = sparse(L_new);
        A_new = sparse(A_new);

        K_new  = d_nbr * W_new;
        y0_new = B_new;

        B_red = SNevo(kappa_new, c, K_new, y0_new);

        recA_conf_i(j+1) = B_red' * L_new * B_red / sum(A_new(:));
        recA_deg_i(j+1)  = sum(A_new(:)) / N_red;
    end

   
    for j = 1:10
        n_rm = 10 * j;

        [A_new, D_new, B_new, kappa_new] = params_filter_B(n_rm, A, D, B_full, kappa);

        N_red = size(D_new,1);
        if N_red == 0
            recB_conf_i(j+1) = NaN;
            recB_deg_i(j+1)  = NaN;
            continue;
        end

        W_new = D_new^(-1) * A_new;
        L_new = D_new - A_new;

        W_new = sparse(W_new);
        L_new = sparse(L_new);
        A_new = sparse(A_new);

        K_new  = d_nbr * W_new;
        y0_new = B_new;

        B_red = SNevo(kappa_new, c, K_new, y0_new);

        recB_conf_i(j+1) = B_red' * L_new * B_red / sum(A_new(:));
        recB_deg_i(j+1)  = sum(A_new(:)) / N_red;
    end

    recA_conf(i,:) = recA_conf_i;
    recA_deg(i,:)  = recA_deg_i;
    recB_conf(i,:) = recB_conf_i;
    recB_deg(i,:)  = recB_deg_i;
end
toc;

save('./result/Compare_ci.mat', ...
     'recA_conf', 'recA_deg', 'recB_conf', 'recB_deg', ...
     'n_rep', 'N', 'alpha');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot: CI vs removed percentage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
load('./result/Compare_ci.mat');
load('color.mat');

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
ax = axes('Parent', gcf);
ax.Position = [0.2, 0.25, 0.6, 0.6];

tick_size  = 18;
label_size = 18;
hold(ax, 'on');

x_vals = 0:10:100;     % number of removed nodes

lA = plot(ax, x_vals, mean(recA_conf, 1, 'omitnan'), ...
    'LineWidth', 4, 'Color', color12, 'DisplayName', 'By Order');
lB = plot(ax, x_vals, mean(recB_conf, 1, 'omitnan'), ...
    'LineWidth', 4, 'Color', color22, 'DisplayName', 'By Disag.');

grid(ax, 'on'); box(ax, 'on');

xlim(ax, [0, 100]);
ylim(ax, [0.15, 0.35]);

xticks(ax, 0:25:100);
xticklabels(ax, {'','1%','2%','3%','4%'});  % 25/2500 = 1%

yticks(ax, 0.2:0.1:0.3);

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);

xlabel('Percentage', 'FontName', 'Arial', 'FontSize', label_size);
% ylabel('Conflict Index', 'FontName', 'Arial', 'FontSize', label_size);

legend([lA, lB], 'FontName', 'Arial', 'FontSize', tick_size, ...
       'Location', 'SouthWest');

title('$\alpha = 5/6$', 'FontSize', label_size, ...
      'FontName', 'Arial', 'Interpreter', 'latex');

set(ax, 'Color', 'none');

saveas(gcf, './fig/fig4_c.svg');

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


function [A_out, D_out, B_out, kappa_out] = params_filter_A(n_rm, A, D, B, kappa)


N     = size(A,1);
A_new = A(n_rm+1:N, n_rm+1:N);
d_vec = A_new * ones(N - n_rm, 1);
D_new = diag(d_vec);

B_new     = B(n_rm+1:N);
kappa_new = kappa(n_rm+1:N);

zero_idx = find(d_vec == 0);
if ~isempty(zero_idx)
    D_new(zero_idx, :) = [];
    D_new(:, zero_idx) = [];
    A_new(zero_idx, :) = [];
    A_new(:, zero_idx) = [];
    B_new(zero_idx)    = [];
    kappa_new(zero_idx)= [];
end

A_out     = A_new;
D_out     = D_new;
B_out     = B_new;
kappa_out = kappa_new;
end


function [A_out, D_out, B_out, kappa_out] = params_filter_B(n_rm, A, D, B, kappa)

N = size(A,1);

loc_conf = (A .* (B' - B).^2) * ones(N,1);
[~, mask] = maxk(loc_conf, n_rm);

A_new = full(A);
A_new(mask, :) = [];
A_new(:, mask) = [];

d_vec = A_new * ones(N - n_rm, 1);
D_new = diag(d_vec);

B_new     = B;
kappa_new = kappa;
B_new(mask)     = [];
kappa_new(mask) = [];

zero_idx = find(d_vec == 0);
if ~isempty(zero_idx)
    D_new(zero_idx, :) = [];
    D_new(:, zero_idx) = [];
    A_new(zero_idx, :) = [];
    A_new(:, zero_idx) = [];
    B_new(zero_idx)    = [];
    kappa_new(zero_idx)= [];
end

A_out     = A_new;
D_out     = D_new;
B_out     = B_new;
kappa_out = kappa_new;
end
