clear; 
close all; 
clc;

year_of_infection = [0,1,2,3,4,5,9,15,16,18,19,25,26,27,28,29,30,31,32,33,34,35,36,37,38]';
Estimated_HIV_incident_numbers = [20.77,65.39,65.13,131.07,131.37,85.37,85.15,49.53,49.89,58.6,59.19,56.40,56.49,53.70,42.21,41.10,40.28,38.89,38.07,38.66,38.95,38.68,37.85,37.30,36.19]';

years_diagnose = [27,28,29,30,31,32,33,34,35,36,37,38]';
HIV_diagnose_cases = [46703,44324,42775,40945,40183,38931,39653,39543,39163,38048,37132,36350]'/1000;

years_AIDS_classification = (4:38)';
AIDS_classification_cases = [11808,18806,27929,34432,41588,47387,57222,69709,73971,67469,64353,56500,45946,39343,37728,38049,36622,36429,36992,35900,33956,32616,31821,31234,30063,27328,25556,24608,23574,19244,18539,18284,17727,17106,16410]'/1000;

d_tau  = 0.52;
kappa  = 52;
dt = d_tau/kappa;
Time_Epidemic  = 39;        % 1981–2020
Treatment_Age  = 20*52;

tforward = 0:dt:Time_Epidemic;
treatment_timeforward = 0:d_tau:Treatment_Age;

base_year = 1981;
real_tforward          = tforward + base_year;
real_tforward_trimmed  = tforward(1:end-1) + base_year;
real_year_of_infection = year_of_infection + base_year;
real_years_diagnose    = years_diagnose    + base_year;
real_years_AIDS_classification = years_AIDS_classification + base_year;

% params = [Lambda, beta, gamma, rho, betad_tilde]
params        = [8.3866, 0.3063, 0.1210, 0.2000, 0.0012];
fitted_params = [33.84, 0.0087, 0.00036, 0.085, 0.45, 0.96];

[HIV_Cases_avg, HIV_Diagnose_avg, AIDS_Cases_avg, ~] = ...
    FiniteDifference_HIV_Between_Host_Model(params, fitted_params, d_tau, kappa, Time_Epidemic, Treatment_Age);

% Use LaTeX so legend math 
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%  STYLE 
fsTick  = 13;     % tick label font
fsLabel = 16;     % axis label font
fsTitle = 19;     % title font
lwAxes  = 3;    % axis line width



lwThin  = 1.4;    % spaghetti lines
lwFit   = 3.0;    % fitted line
msCDC   = 7;      % CDC marker size

cLines  = [0.62 0.36 0.86];  % purple (individuals)
cFit    = [0.00 0.35 0.80];  % deep blue (fit)
cCDC    = [0.80 0.12 0.12];  % dark red (CDC)

%  LOAD INDIVIDUAL PARAM SETS 
filename = 'Parameters_pop.xlsx';
T = readtable(filename, 'ReadVariableNames', true, 'PreserveVariableNames', true);

lambda_vals  = getVar(T, {'lambda'});
k_vals       = getVar(T, {'beta','k'});
delta_1_vals = getVar(T, {'d','delta_1'});
delta_2_vals = getVar(T, {'delta','delta_2'});
p_vals       = getVar(T, {'pi','p'});
c_vals       = getVar(T, {'c'});

num_individuals = numel(lambda_vals);

%  RUN INDIVIDUALS 
Nt = numel(tforward);
Incidence  = nan(Nt-1, num_individuals);
Diagnoses  = nan(Nt,   num_individuals);
AIDSflows  = nan(Nt-1, num_individuals);

for i = 1:num_individuals
    wh = [lambda_vals(i), k_vals(i), delta_1_vals(i), ...
          delta_2_vals(i), p_vals(i), c_vals(i)];
    [HIV_Cases_i, HIV_Dx_i, AIDS_Cases_i, ~] = ...
        FiniteDifference_HIV_Between_Host_Model(params, wh, d_tau, kappa, Time_Epidemic, Treatment_Age);
    Incidence(:,i) = HIV_Cases_i;
    Diagnoses(:,i) = HIV_Dx_i;
    AIDSflows(:,i) = AIDS_Cases_i;
end


% Incidence 
figure('Position',[60,80,1100,460]); hold on; box on; grid on;
plot(real_tforward_trimmed, Incidence, '-', 'Color', cLines, 'LineWidth', lwThin);
plot(real_tforward_trimmed, HIV_Cases_avg, '-', 'Color', cFit, 'LineWidth', lwFit);
plot(real_year_of_infection, Estimated_HIV_incident_numbers, 'o', ...
    'MarkerEdgeColor','k','MarkerFaceColor',cCDC,'MarkerSize',msCDC);
xlim([1981 2020]); xticks(1985:5:2020);

% Bold title/labels even with LaTeX interpreter
title('\textbf{Estimated Annual HIV Incidence}');
xlabel('\textbf{Year}');
ylabel('\textbf{HIV incidence (in thousands)}');




ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
beautify(ax, fsTick, fsLabel, fsTitle, lwAxes);
ax.TickLabelInterpreter = 'tex';   % allow tick labels to be bold
ax.FontWeight = 'bold';
ax.LineWidth = 3.5;      % Thicker frame
ax.TickLength = [0.025 0.02];  % Longer ticks

% Bold legend (math-less)
hI1 = plot(nan,nan,'-','Color',cLines,'LineWidth',lwThin);
hI2 = plot(nan,nan,'-','Color',cFit,  'LineWidth',lwFit);
hI3 = plot(nan,nan,'o','MarkerEdgeColor','k','MarkerFaceColor',cCDC,'MarkerSize',msCDC);
leg = legend([hI1,hI2,hI3], ...
    {'\textbf{Individuals}','\textbf{Fitted prediction}','\textbf{CDC data}'}, ...
    'Location','northwest','Interpreter','latex');
set(leg,'Box','off','FontSize',12);

% Diagnoses 
figure('Position',[80,110,1100,460]); hold on; box on; grid on;
plot(real_tforward, Diagnoses, '-', 'Color', cLines, 'LineWidth', lwThin);
plot(real_tforward, HIV_Diagnose_avg, '-', 'Color', cFit, 'LineWidth', lwFit);
plot(real_years_diagnose, HIV_diagnose_cases, 'o','MarkerEdgeColor','k','MarkerFaceColor',cCDC,'MarkerSize',msCDC);
xlim([1981 2020]); xticks(1985:5:2020);

title('\textbf{Annual HIV Diagnoses}');
xlabel('\textbf{Year}');
ylabel('\textbf{HIV diagnoses (in thousands)}');

ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
beautify(ax, fsTick, fsLabel, fsTitle, lwAxes);
ax.TickLabelInterpreter = 'tex';  % bold ticks
ax.FontWeight = 'bold';
ax.LineWidth = 3.5;      % Thicker frame
ax.TickLength = [0.025 0.02];  % Longer ticks

hD1 = plot(nan,nan,'-','Color',cLines,'LineWidth',lwThin);
hD2 = plot(nan,nan,'-','Color',cFit,  'LineWidth',lwFit);
hD3 = plot(nan,nan,'o','MarkerEdgeColor','k','MarkerFaceColor',cCDC,'MarkerSize',msCDC);
leg = legend([hD1,hD2,hD3], ...
    {'\textbf{Individuals}','\textbf{Fitted prediction}','\textbf{CDC data}'}, ...
    'Location','northwest','Interpreter','latex');
set(leg,'Box','off','FontSize',12);

%  AIDS classifications 
figure('Position',[100,140,1100,460]); hold on; box on; grid on;
plot(real_tforward_trimmed, AIDSflows, '-', 'Color', cLines, 'LineWidth', lwThin);
plot(real_tforward_trimmed, AIDS_Cases_avg, '-', 'Color', cFit, 'LineWidth', lwFit);
plot(real_years_AIDS_classification, AIDS_classification_cases, 'o', ...
    'MarkerEdgeColor','k','MarkerFaceColor',cCDC,'MarkerSize',msCDC);
xlim([1981 2020]); xticks(1985:5:2020);

title('\textbf{Annual AIDS Classifications}');
xlabel('\textbf{Year}');
ylabel('\textbf{AIDS classifications (in thousands)}');

ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
beautify(ax, fsTick, fsLabel, fsTitle, lwAxes);
ax.TickLabelInterpreter = 'tex';  % bold ticks
ax.FontWeight = 'bold';
ax.LineWidth = 3.5;      % Thicker frame
ax.TickLength = [0.025 0.02];  % Longer ticks

hA1 = plot(nan,nan,'-','Color',cLines,'LineWidth',lwThin);
hA2 = plot(nan,nan,'-','Color',cFit,  'LineWidth',lwFit);
hA3 = plot(nan,nan,'o','MarkerEdgeColor','k','MarkerFaceColor',cCDC,'MarkerSize',msCDC);
leg = legend([hA1,hA2,hA3], ...
    {'\textbf{Individuals}','\textbf{Fitted prediction}','\textbf{CDC data}'}, ...
    'Location','northwest','Interpreter','latex');
set(leg,'Box','off','FontSize',12);



odeopts = odeset('RelTol',1e-8,'AbsTol',1e-12);
init    = [167; 85; 69643];           % [T0; Ti0; V0]
x_weeks = treatment_timeforward;      % τ in weeks

beta_d_all = nan(numel(x_weeks), num_individuals);

for ii = 1:num_individuals
    wh = [lambda_vals(ii), k_vals(ii), delta_1_vals(ii), ...
          delta_2_vals(ii), p_vals(ii), c_vals(ii)];
    [~, y_sol] = ode23s(@(t,y) HIV_Within_Host(y,wh), x_weeks, init, odeopts);
    V = y_sol(:,3);
    beta_week = params(5) * V;        
    beta_d_all(:,ii) = 52 * beta_week; % plot in PER YEAR
end


%  Plot β_d(τ) 
figure('Position',[120,170,900,420]); 
hold on; box on; grid on;

% Plot individual within-host trajectories only
plot(x_weeks, beta_d_all, '-', 'Color', [0.62 0.36 0.86], 'LineWidth', 0.6);

% Axis formatting (match other figures)
title('\textbf{Within-host Transmission Rate $\beta_d(\tau)$}','Interpreter','latex');
xlabel('\textbf{Treatment age (weeks)}','Interpreter','latex');
ylabel('\textbf{Transmission rate $\beta_d$ }','Interpreter','latex');

xlim([0, Treatment_Age]);        % full range 0 → 1040 weeks
ax = gca;

%  Match style of other figures ---
beautify(ax, fsTick, fsLabel, fsTitle, lwAxes);
ax.TickLabelInterpreter = 'tex';
ax.FontWeight = 'bold';
ax.LineWidth = 2;      % Thicker frame
ax.TickLength = [0.025 0.02];  % Longer ticks

% Grid: only horizontal lines
ax.XGrid = 'off';
ax.YGrid = 'on';


ax.Box = 'on';       
ax.Layer = 'top';   



function [HIV_Cases, HIV_Diagnose, AIDS_Cases, AIDS_deaths] = ...
    FiniteDifference_HIV_Between_Host_Model(params, within_params, d_tau, kappa, Time_Epidemic, Treatment_Age)

    dt = d_tau/kappa;
    tforward = 0:dt:Time_Epidemic;
    treatment_timeforward = 0:d_tau:Treatment_Age;

    Lambda = params(1);
    mu     = 0.02;
    beta   = params(2);
    gamma  = params(3);
    rho    = params(4);
    alpha  = 0.002;

    % Within-host
    T0=167; T_i0=85; V0=69643; y0=[T0; T_i0; V0];
    odeopts = odeset('RelTol',1e-8,'AbsTol',1e-12);
    [~, y_sol] = ode23s(@(t,y) HIV_Within_Host(y, within_params), treatment_timeforward, y0, odeopts);
    beta_d_vec = params(5) * y_sol(:,3);   % linear form inside the FD model

    % Finite difference state variables
    N_steps = numel(tforward);
    M = numel(treatment_timeforward);

    S = zeros(N_steps,1); I = zeros(N_steps,1); A = zeros(N_steps,1);
    d_old = zeros(M,1);   d_new = zeros(M,1);   D = zeros(N_steps,1);

    S(1)=820; I(1)=42; A(1)=0.337;
    D(1)=trapz(d_old)*d_tau;

    Bd_Integral = zeros(N_steps,1);
    N_total     = zeros(N_steps,1);
    N_total(1)  = S(1)+I(1)+D(1);

    for k=1:N_steps-1
        N_total(k) = S(k)+I(k)+D(k);
        Den = max(N_total(k), 1e-12);
        Bd_Integral(k) = d_tau*trapz(beta_d_vec.*d_old)/Den;

        S(k+1) = (S(k)+dt*Lambda) / (1+dt*(beta*I(k)/Den+mu));
        I(k+1) = (I(k)+dt*(beta*S(k+1)*I(k)/Den + S(k+1)*Bd_Integral(k))) / (1+(mu+alpha+rho)*dt);

        d_new(1)    = (1/kappa)*rho*I(k+1);
        d_new(2:end)= d_old(1:end-1) ./ (1+dt*mu);

        A(k+1) = (A(k)+dt*(gamma*I(k+1))) / (1+dt*(mu+alpha));
        D(k+1) = trapz(d_new)*d_tau;

        d_old = d_new;
    end

    HIV_Cases    = beta*S(1:end-1).*I(1:end-1)./max(N_total(1:end-1),1e-12) + S(1:end-1).*Bd_Integral(1:end-1);
    HIV_Diagnose = rho*I;
    AIDS_Cases   = gamma*I(1:end-1);
    AIDS_deaths  = (mu + alpha)*A;
end

function dy = HIV_Within_Host(y, p)
    T=y(1); Ti=y(2); V=y(3);
    lambda=p(1); k=p(2); delta_1=p(3); delta_2=p(4); P=p(5); c=p(6);
    T = max(T,0); Ti = max(Ti,0); V = max(V,0);  % positivity guard (optional)
    dy=zeros(3,1);
    dy(1)= lambda - k*T*V - delta_1*T;
    dy(2)= k*T*V - delta_2*Ti;
    dy(3)= P*Ti   - c*V;
end

function v = getVar(T, aliases)
    vnames = string(T.Properties.VariableNames);
    lowerV = lower(vnames);
    found = [];
    for a = string(aliases)
        j = find(lowerV == lower(a), 1);
        if ~isempty(j), found = j; break; end
    end
    if isempty(found)
        error('Missing required column. Tried aliases: %s', strjoin(aliases, ', '));
    end
    v = T.(vnames(found));
end

function beautify(ax, fsTick, fsLabel, fsTitle, lwAxes)
    if nargin < 1 || isempty(ax), ax = gca; end
    if nargin < 2 || isempty(fsTick),  fsTick  = 13; end
    if nargin < 3 || isempty(fsLabel), fsLabel = 16; end
    if nargin < 4 || isempty(fsTitle), fsTitle = 19; end
    if nargin < 5 || isempty(lwAxes),  lwAxes  = 1.8; end

    set(ax, ...
        'FontSize', fsTick, ...
        'LineWidth', lwAxes, ...
        'TickDir', 'out', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.015 0.015], ...
        'Box', 'on', ...
        'GridAlpha', 0.15, 'MinorGridAlpha', 0.08);

    ax.Title.FontSize   = fsTitle;
    ax.Title.FontWeight = 'bold';
    ax.XLabel.FontSize  = fsLabel;
    ax.YLabel.FontSize  = fsLabel;
end
