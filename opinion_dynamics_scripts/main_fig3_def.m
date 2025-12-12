%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig4_networks.m
%
% This script generates Fig.3(d)-(f) for article:
%
% Opinion polarization and its connected disagreement: Modeling and modulation
%
% Author: Xuzhe Qian
% Date: 2025/12/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------------------------------------------------
%  Common settings
% -------------------------------------------------------------------------
clear; clc; close all;

%% ------------------------------------------------------------------------
%  WS network: empirical eigenvalue distribution
% -------------------------------------------------------------------------

N       = 2500;
ad      = 4;
p       = 0.1;
n_total = 100;

rec_sp  = zeros(n_total, N);
rec_min = zeros(n_total, 1);

tic;
parfor n = 1:n_total

    A   = sw_net('N', N, 'k', ad, 'p', p);
    A   = sparse(A);
    deg = A * ones(N,1);

    W     = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L_bar = eye(N) - W;

    sp = eig(L_bar);
    rec_sp(n,:)  = sp;
    rec_min(n)   = min(real(sp(2:end)));
end
toc;

save('./result/sp_ws.mat', 'rec_sp', 'rec_min', 'n_total', 'N');
rec_sp = rec_sp(:);

load('./result/sp_ws.mat');  % rec_sp, rec_min, n_total, N

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
ax = axes('Parent', gcf); ax.Position = [0.2, 0.25, 0.6, 0.6];
tick_size  = 18;
label_size = 18;
hold(ax, 'on');

load('color.mat');
rec_sp = real(rec_sp);

h = histogram(ax, rec_sp, 0:0.01:2);
h.FaceColor = color1;
h.EdgeColor = 'none';

eig_min = mean(rec_min);
plot(ax, [eig_min, eig_min], [0, 37.5 * n_total], ...
    'LineWidth', 4, 'LineStyle', ':', 'Color', color12);

grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, 1]);
ylim(ax, [0, 37.5 * n_total]);

xticks(ax, 0:0.25:1);
yticks(ax, 0:12.5*n_total:37.5*n_total);
yticklabels(ax, {'0','0.5','1','1.5'});

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);

xlabel('Eigenvalue {\lambda}', 'FontName', 'Arial', 'FontSize', label_size);
ylabel('Empirical Density',   'FontName', 'Arial', 'FontSize', label_size);

saveas(gcf, './fig/fig3_ws_sp.svg');


%% ------------------------------------------------------------------------
%  PS network: empirical eigenvalue distribution
% -------------------------------------------------------------------------

clear; clc; close all;

n_total   = 100;
N         = 2500;
m         = 2;
T         = 0.2;
gamma     = 2.1;
distr     = 0;
plot_flag = 0;

rec_sp  = zeros(n_total, N);
rec_min = zeros(n_total, 1);

tic;
parfor n = 1:n_total
    [A, coords, comm, d] = nPSO_model(N, m, T, gamma, distr, plot_flag); 
    A   = sparse(A);
    deg = A * ones(N,1);

    W     = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L_bar = eye(N) - W;

    sp = eig(L_bar);
    rec_sp(n,:) = real(sp);
    rec_min(n)  = min(real(sp(2:end)));
end
toc;

rec_sp = real(rec_sp(:));
save('./result/sp_ps.mat', 'rec_sp', 'rec_min', 'n_total', 'N');

%% Plot: PS network spectrum
load('./result/sp_ps.mat');

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
ax = axes('Parent', gcf); ax.Position = [0.2, 0.25, 0.6, 0.6];
tick_size  = 18;
label_size = 18;
hold(ax, 'on');

load('color.mat');

h = histogram(ax, real(rec_sp), 0:0.01:2);
h.FaceColor = color1;
h.EdgeColor = 'none';

eig_min = mean(rec_min);
plot(ax, [eig_min, eig_min], [0, 37.5 * n_total], ...
    'LineWidth', 4, 'LineStyle', ':', 'Color', color12);

grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, 1]);
ylim(ax, [0, 37.5 * n_total]);

xticks(ax, 0:0.25:1);
yticks(ax, 0:12.5*n_total:37.5*n_total);
yticklabels(ax, {'0','0.5','1','1.5'});

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);

xlabel('Eigenvalue {\lambda}', 'FontName', 'Arial', 'FontSize', label_size);
ylabel('Empirical Density',   'FontName', 'Arial', 'FontSize', label_size);

saveas(gcf, './fig/fig3_ps_sp.svg');


%% ------------------------------------------------------------------------
%  BA network: empirical eigenvalue distribution
% -------------------------------------------------------------------------

clear; clc; close all;

N       = 2500;
ad      = 4;         
n_total = 100;

rec_sp  = zeros(n_total, N);
rec_min = zeros(n_total, 1);

tic;
parfor n = 1:n_total
    A   = ba_net('N', N, 'm0', 2, 'm', 2);
    A   = sparse(A);
    deg = A * ones(N,1);
    W     = sparse(spdiags(deg.^(-1), 0, N, N)) * A;
    L_bar = eye(N) - W;
    sp = eig(L_bar);
    rec_sp(n,:) = real(sp);
    rec_min(n)  = min(real(sp(2:end)));
end
toc;

rec_sp = real(rec_sp(:));
save('./result/sp_ba.mat', 'rec_sp', 'rec_min', 'n_total', 'N');

%% Plot: BA network spectrum
load('./result/sp_ba.mat');

figure('OuterPosition',[800, 200, 420, 350], 'visible', 'on');
ax = axes('Parent', gcf); ax.Position = [0.2, 0.25, 0.6, 0.6];
tick_size  = 18;
label_size = 18;
hold(ax, 'on');

load('color.mat');

h = histogram(ax, real(rec_sp), 0:0.01:2);
h.FaceColor = color1;
h.EdgeColor = 'none';

eig_min = mean(rec_min);
plot(ax, [eig_min, eig_min], [0, 37.5 * n_total], ...
    'LineWidth', 4, 'LineStyle', ':', 'Color', color12);

grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, 1]);
ylim(ax, [0, 37.5 * n_total]);

xticks(ax, 0:0.25:1);
yticks(ax, 0:12.5*n_total:37.5*n_total);
yticklabels(ax, {'0','0.5','1','1.5'});

set(ax.XAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);
set(ax.YAxis, 'FontSize', tick_size, 'FontName', 'Arial', 'LineWidth', 4/3);

xlabel('Eigenvalue {\lambda}', 'FontName', 'Arial', 'FontSize', label_size);
ylabel('Empirical Density',   'FontName', 'Arial', 'FontSize', label_size);

saveas(gcf, './fig/fig3_ba_sp.svg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
