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

clear; close all;

%% ------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------
kappa = 1;
K = 0;
y0 = 1;                     % initial opinion
kappa  = 3;                 % sensitivity parameter
F_step = 0.005;             % step size of external forcing scan
tspan = 0:50:1000;          % ODE integration interval
ode_opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% ------------------------------------------------------------------------
%  Pre-allocation helper
% -------------------------------------------------------------------------
scan_range = 200;                     % corresponds to 0 → ±1 in F_ext
xN = -1:F_step:1;

%% ------------------------------------------------------------------------
%  Compute curves for alpha [0, 2/3]
% -------------------------------------------------------------------------
alphas = [0, 2/3];               
rec_curves = cell(length(alphas),1);

for idx = 1:length(alphas)
    alpha = alphas(idx);
    c = alpha;
    d = (1-alpha);

    rec = zeros(1,scan_range);
    for i = 1:scan_range
        F_ext = d * (i*F_step);
        [~,Y] = ode45(@(t,y) SN(t,y,kappa,c,K,F_ext), tspan, y0, ode_opts);
        rec(i) = Y(end,1);
    end
    rec_curves{idx} = [-rec(end:-1:1), 0, rec];   % symmetric extension
end

rec15 = rec_curves{1};
rec35 = rec_curves{2};

%% ------------------------------------------------------------------------
%  Compute curves for alpha [5/6]
% -------------------------------------------------------------------------
alpha = 5/6;
c = alpha;
d = (1-alpha);

rec45 = zeros(1,2*scan_range+1);
cnt = 0;
for i = -scan_range:scan_range
    cnt = cnt + 1;
    F_ext = d * (i*F_step);
    [~,Y] = ode45(@(t,y) SN(t,y,kappa,c,K,F_ext), tspan, y0, ode_opts);
    rec45(cnt) = Y(end,1);
end

d45 = -1:F_step:1;
transition_index = sum(rec45 < 0);
s_val = d45(transition_index+1);

% Stable Branches
rec45_p = rec45(transition_index+1:end);
rec45_n = -rec45(end:-1:transition_index+1);

d45_p = d45(transition_index+1:end);
d45_n = -d45(end:-1:transition_index+1);

% Unstable Branch
b_vals = s_val : F_step : -s_val;  
rec45_u = zeros(size(b_vals));
for j = 1:length(b_vals)
    b = b_vals(j);
    x_guess = 0;
    fun = @(x) x - 2/(1 + exp(kappa*(-c*x + b*d))) + 1;
    [x_sol, ~] = fsolve(fun, x_guess, optimset('Display','off'));
    rec45_u(j) = -x_sol;
end
d45_u = b_vals;

%% ------------------------------------------------------------------------
%  Plotting
% -------------------------------------------------------------------------
load('color.mat');  

fig_width = 460;
fig_height = 360;
label_size = 18;
tick_size  = 18;

figure('OuterPosition',[300,300,fig_width,fig_height]);
hold on;

h1 = plot(xN, rec15, 'LineWidth', 4, 'Color', color12, ...
          'DisplayName', '\alpha = 0');
h2 = plot(xN, rec35, 'LineWidth', 4, 'Color', color12, ...
          'DisplayName', '\alpha = 2/3');

h3p = plot(d45_p, rec45_p, 'LineWidth', 4, 'Color', color12, ...
           'DisplayName', '\alpha = 5/6');
h3n = plot(d45_n, rec45_n, 'LineWidth', 4, 'Color', color12);
h3u = plot(d45_u, rec45_u, ':', 'LineWidth', 4, 'Color', color12);

h1.Color(4)  = 0.4;
h2.Color(4)  = 0.7;
h3p.Color(4) = 1;
h3n.Color(4) = 1;
h3u.Color(4) = 1;

xlim([-0.75, 0.75]);
xticks(-0.5:0.5:0.5);
ylim([-1, 1]);
yticks(-1:0.5:1);
yticklabels({'-1','','0','','1'});

ax = gca;
set(ax.XAxis,'FontSize',tick_size,'FontName','Arial','LineWidth',4/3);
set(ax.YAxis,'FontSize',tick_size,'FontName','Arial','LineWidth',4/3);
set(gca,'FontSize',label_size,'Color','none','LineWidth',4/3);

grid on; box on;
legend([h1,h2,h3p],'FontSize',18,'FontName','Arial','Location','East');

xlabel('Average neighboring opinion','FontSize',label_size,'FontName','Arial');
ylabel('Opinion','FontSize',label_size,'FontName','Arial');

saveas(gcf,'./fig/fig1_b.svg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = SN(~, y, kappa, c, K, F_ext)

B = y;
E = c*B + K*B + F_ext;

dydt = -B + 2./(1 + exp(-kappa.*E)) - 1;
end

