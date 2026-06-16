clear all; 
close all; 
clc;

% Fitted parameters we found  
Lambda = 8.38663746616526;
beta = 0.306293687587148;
gamma = 0.120967014736454;
rho = 0.200000005673198;  
beta_d0 = 0.00119503317441228; % fitted tilde beta_d up to 2019
params = [Lambda, beta, gamma, rho, beta_d0];

% Simulation time setup
base_year = 1981;
dt = 0.01; 
Time_Epidemic = 49; % extend to 2030
transition_start = 2019;
transition_end = 2030;

% Target goals (in thousands)
goal2025 = 9.3;
goal2030 = 3.0;
goal2025_diag = 9.588;
goal2030_diag = 3.0;

% CDC Data for plotting
year_of_infection = [0, 1, 2, 3, 4, 5, 9, 15, 16, 18, 19, 25, 26, 27, 28, ...   
  29, 30, 31, 32, 33, 34, 35, 36, 37, 38]';

Estimated_HIV_incident_numbers = [20.77, 65.39, 65.13, 131.07, 131.37,...
    85.37, 85.15, 49.53, 49.89, 58.6, 59.19, 56.40, 56.49, 53.70,...
    42.21, 41.10, 40.28, 38.89, 38.07, 38.66, 38.95,...
    38.68, 37.85, 37.30, 36.19]';

year_of_infection_post = [39, 40, 41]';   % 2020–2022

Estimated_HIV_incident_numbers_post = [34.2, 32.7, 31.8]'; % in 1000s

years_diagnose = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]';
HIV_diagnose_cases = [46703, 44324, 42775, 40945, 40183, 38931, 39653,...
39543, 39163, 38048, 37132, 36350]'/1000;

years_diagnose_post = [39, 40, 41, 42, 43]'; % 2020–2024

HIV_diagnose_cases_post = [30.403, 35.763, 37.660, 38.793, 34.788]'; % in 1000s

years_AIDS_classification = (4:38)';

AIDS_classidication_cases = [11808, 18806, 27929, 34432, 41588, 47387, 57222,...
69709, 73971, 67469, 64353, 56500, 45946, 39343, 37728, 38049, 36622,...
36429, 36992, 35900, 33956, 32616, 31821, 31234, 30063, 27328, 25556,...
24608, 23574, 19244, 18539, 18284, 17727, 17106, 16410]'/1000;

years_AIDS_death = (4:38)';

AIDS_death_cases = [6712, 11727, 15703, 19980, 26438, 30051, 35145, 39861,...
44158, 49396, 50468, 37450, 22056, 18857, 18587, 18275, 18397, 18339,...
18202, 17799, 17470, 17094, 16288, 15361, 14847, 13703, 13330, 13041,...
12977, 13060, 12779, 12945, 12825, 12649, 12471]'/1000;

% Objective function for post-2019 beta_d_tilde
objfun = @(beta_d_post) ...
    10*(get_incidence_at_year(2025, params, dt, Time_Epidemic, base_year, beta_d0, beta_d_post, transition_start, transition_end) - goal2025)^2 + ...
    20*(get_incidence_at_year(2030, params, dt, Time_Epidemic, base_year, beta_d0, beta_d_post, transition_start, transition_end) - goal2030)^2 + ...
    10*(get_diagnosis_at_year(2025, params, dt, Time_Epidemic, base_year, beta_d0, beta_d_post, transition_start, transition_end) - goal2025_diag)^2 + ...
    20*(get_diagnosis_at_year(2030, params, dt, Time_Epidemic, base_year, beta_d0, beta_d_post, transition_start, transition_end) - goal2030_diag)^2;

% Optimize only beta_d_post 
nonlcon = @(x) beta_constraints(x, beta);  % beta is a known parameter
beta_d_post_opt = fmincon(objfun, beta_d0, [], [], [], [], 0, 10*beta_d0, nonlcon);

fprintf('Optimized post-2019 beta_d_tilde = %.6f\n', beta_d_post_opt);

% R0: constant <=2019; varies with \tilde{beta}_d(t) after 2019 
[R0_years, R0_series, R0_base, R0_coef, R0_pre2019] = ...
    compute_R0_series(params, dt, Time_Epidemic, base_year, ...
                      beta_d0, beta_d_post_opt, transition_start, transition_end);

% (Optional) print a few reference values
R0_2018  = R0_series(find(abs(R0_years-2018)==min(abs(R0_years-2018)),1));
R0_2025  = R0_series(find(abs(R0_years-2025)==min(abs(R0_years-2025)),1));
R0_2030  = R0_series(find(abs(R0_years-2030)==min(abs(R0_years-2030)),1));
fprintf('R0 base term (beta/(mu+gamma+rho))         = %.6f\n', R0_base);
fprintf('R0 diagnosed-term coefficient              = %.6f\n', R0_coef);
fprintf('R0 (<=2019, constant total)                = %.6f\n', R0_pre2019);
fprintf('R0(2018) = %.6f,  R0(2025) = %.6f,  R0(2030) = %.6f\n', R0_2018, R0_2025, R0_2030);




% Generate full projection w/ optimized beta_d for plotting
[years, HIV_Inc, HIV_Diag] = project_full(params, dt, Time_Epidemic, base_year, beta_d0, beta_d_post_opt, transition_start, transition_end);

% Time vectors for plotting
tforward = 0:dt:(2030-base_year);
real_t_inc = tforward(1:end-1) + base_year; 
real_t_full = tforward + base_year;



% Plot: HIV Incidence 
figure('Position', [100, 100, 1000, 600]); hold on;

xline(2019, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off');

scatter(2025, goal2025, 100, [0 0.6 0], 'v', 'filled', 'HandleVisibility', 'off');

scatter(2030, goal2030, 100, [0 0 0.8], 'v', 'filled', 'HandleVisibility', 'off');

p1 = plot(real_t_inc, HIV_Inc, '-b', 'LineWidth', 2);

p2 = scatter(year_of_infection + base_year, Estimated_HIV_incident_numbers, 80, [0.8 0 0.8], 'filled');

real_year_inf_post = year_of_infection_post + base_year;

p3 = scatter(real_year_inf_post, Estimated_HIV_incident_numbers_post, ...
             80, [0.85 0 0], 'filled');

legend([p1 p2 p3], ...
   {'Model predicted HIV incidences', ...
    'CDC reported incidences (≤2019)', ...
    'CDC reported incidences (2020+)'}, ...
   'Location', 'northwest', 'AutoUpdate','off');

xlabel('Years of HIV Epidemic in USA', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('HIV Cases Per Year (in 1000s)', 'FontWeight', 'bold', 'FontSize', 12);
title('HIV Incidences Projected to 2030');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
box on; grid on; xlim([1981, 2032]); 
yl = ylim;
text(2013 + 0.2, yl(2) * 0.75, ...
'Ending the HIV Epidemic Initiated', ...
'FontSize', 14, 'Color', [0 0.5 0]);
% Goal annotations
arrow_offset = 40;
ax = gca;
ax_pos = get(ax, 'Position');
xlims = xlim(ax);
ylims = ylim(ax);
x2fig = @(x) ax_pos(1) + (x - xlims(1)) / diff(xlims) * ax_pos(3);
y2fig = @(y) ax_pos(2) + (y - ylims(1)) / diff(ylims) * ax_pos(4);


annotation('textarrow', ...
    [x2fig(2025 - 2.5), x2fig(2025)], ...
    [y2fig(goal2025 + arrow_offset), y2fig(goal2025)], ...
    'String', '2025 Year Target: 9.3k', ...
    'FontSize', 12, 'Color', [0 0.5 0]);
annotation('textarrow', ...
    [x2fig(2030 - 3.6), x2fig(2030)], ...
    [y2fig(goal2030 + arrow_offset), y2fig(goal2030)], ...
    'String', '2030 Year Target: 3k', ...
    'FontSize', 12, 'Color', [0 0 0.6]);

% Plot: HIV Diagnoses 
figure('Position', [100, 100, 1000, 500]); 

hold on;

p1 = plot(real_t_inc, HIV_Diag, '-b', 'LineWidth', 2);

p2 = scatter(years_diagnose + base_year, HIV_diagnose_cases, 80, [0.8 0 0.8], 'filled');
xline(2019, '--k', 'LineWidth', 1.5);

real_year_diag_post = years_diagnose_post + base_year;
p3 = scatter(real_year_diag_post, HIV_diagnose_cases_post, ...
              80, 'o', 'MarkerFaceColor',[0.85 0 0], 'MarkerEdgeColor','k');

scatter(2025, goal2025_diag, 100, [0 0.6 0], 'v', 'filled', 'HandleVisibility', 'off');
scatter(2030, goal2030_diag, 100, [0 0 0.8], 'v', 'filled', 'HandleVisibility', 'off');

legend([p1 p2 p3], ...
   {'Model predicted HIV diagnoses', ...
    'CDC reported diagnoses (≤2019)', ...
    'CDC reported diagnoses (2020+)'}, ...
   'Location','northwest','AutoUpdate','off');

xlabel('Years of HIV Epidemic in USA', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('HIV Diagnoses Per Year (in 1000s)', 'FontWeight', 'bold', 'FontSize', 12);

title('HIV Diagnoses Projected to 2030');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
box on; 
grid on;
xlim([1981, 2032]);
yl = ylim;
text(2013 + 0.2, yl(2) * 0.75, ...
'Ending the HIV Epidemic Initiated', ...
'FontSize', 14, 'Color', [0 0.5 0]);
ax = gca;
ax_pos = get(ax, 'Position');
xlim_ax = xlim(ax); ylim_ax = ylim(ax);
x2fig = @(x) ax_pos(1) + (x - xlim_ax(1)) / diff(xlim_ax) * ax_pos(3);
y2fig = @(y) ax_pos(2) + (y - ylim_ax(1)) / diff(ylim_ax) * ax_pos(4);
arrow_offset = 40;

annotation('textarrow', ...
[x2fig(2025 - 2.5), x2fig(2025)], ...
[y2fig(goal2025_diag + arrow_offset), y2fig(goal2025_diag)], ...
'String', '2025 Year Target: 9.6k', ...
'FontSize', 12, 'Color', [0 0.5 0]);

annotation('textarrow', ...
[x2fig(2030 - 3.6), x2fig(2030)], ...
[y2fig(goal2030_diag + arrow_offset), y2fig(goal2030_diag)], ...
'String', '2030 Year Target: 3k', ...
'FontSize', 12, 'Color', [0 0 0.6]);

% Plot: beta_d(t) 
tt = 1981:transition_end;   % transition_end = 2030

beta_d_values = arrayfun(@(t) ...
    (t <= transition_start) * beta_d0 + ...
    (t > transition_start & t <= transition_end) * ...
        ((beta_d_post_opt - beta_d0) / (transition_end - transition_start) ...
        * (t - transition_start) + beta_d0), ...
    tt);

figure('Position', [100, 100, 1200, 600]);
plot(tt, beta_d_values, '-b', 'LineWidth', 2);
% xlabel('Year', 'FontWeight', 'bold', 'FontSize', 14);
% ylabel('$\tilde{\beta}_d(t)$', 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'latex');

% title('Time-Dependent Transmission Rate $\tilde{\beta}_d(t)$', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on; xlim([1981, 2032]); ylim ([0, 0.002]); yticks(0:0.05:0.2);
set(gca, 'FontSize', 13, 'LineWidth', 1.5); box on;


xlabel('Year','FontWeight','bold','FontSize',16);

% ylabel('$\tilde{\beta}_d(t)$', ...
%     'Interpreter','latex', ...
%     'FontSize',18);
% 
% title('Time-Dependent Transmission Rate $\tilde{\beta}_d(t)$', ...
%     'Interpreter','latex', ...
%     'FontSize',18);

ylabel('$\tilde{\beta}_d(t)$', ...
    'Interpreter','latex', ...
    'FontSize',22);

title('Time-Dependent Transmission Rate $\tilde{\beta}_d(t)$', ...
    'Interpreter','latex', ...
    'FontSize',20);


% Figure: R0 over time 
figure('Position',[100,100,950,380]); hold on;

% ensure column vectors for plotting
plot(R0_years(:), R0_series(:), '-b', 'LineWidth', 2);
xline(2019, '--k', 'LineWidth', 1.5);
yline(1, '--', 'LineWidth', 1.2, 'Color', [0 0.6 0]);

xlabel('Year', 'FontWeight','bold', 'FontSize',12);
ylabel('Basic reproduction number  R_0', 'FontWeight','bold', 'FontSize',12);
box on
hR0 = plot(R0_years(:), R0_series(:), '-b', 'LineWidth', 2); hold on;
hCP  = xline(2019, '--', 'LineWidth', 1.2, 'Color', [0 0 0]);
hThr = yline(1,    '--', 'LineWidth', 1.2, 'Color', [0 0.6 0]);

title({'Refitting $\tilde{\beta}_d$ after 2019', ...
       'Basic reproduction number $R_0(t)$'}, ...
      'Interpreter','latex','FontSize',18);

lgd = legend([hR0 hCP hThr], ...
    {'$R_0(t)$', 'Change point: 2019', 'Threshold $R_0=1$'}, ...
    'Location','northeast', 'Interpreter','latex', 'Box','off');
lgd.FontSize = 15;

% (Optional) tweak legend font size



grid on; xlim([1981, 2032]);
ylim([0, max(R0_series)*1.1]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
hold off;



% Helper Functions 
function inc = get_incidence_at_year(year_target, params, dt, T_end, base_year, beta_d0, beta_d_post, transition_start, transition_end)
    [years, HIV_Inc, ~] = project_full(params, dt, T_end, base_year, beta_d0, beta_d_post, transition_start, transition_end);
    [~, idx] = min(abs(years - year_target));
    inc = HIV_Inc(idx);
end

function diag = get_diagnosis_at_year(year_target, params, dt, T_end, base_year, beta_d0, beta_d_post, transition_start, transition_end)
    [years, ~, HIV_Diag] = project_full(params, dt, T_end, base_year, beta_d0, beta_d_post, transition_start, transition_end);
    [~, idx] = min(abs(years - year_target));
    diag = HIV_Diag(idx);
end

function [years, HIV_Inc, HIV_Diag] = project_full(params, dt, T_end, base_year, beta_d0, beta_d_post, t_start_transition, t_end_transition)
    
    
    Lambda = params(1);
    beta = params(2);
    gamma = params(3);
    rho = params(4);

    d_tau = 0.52;
    kappa = 52;
    Treatment_Age = 20*52;

    tforward = 0:dt:T_end;
    treatment_timeforward = 0:d_tau:Treatment_Age;
    
    % Within-host parameters (fixed)
    fitted_params = [33.84, 0.0087, 0.00036, 0.085, 0.45, 0.96];
    lambda_w = fitted_params(1);
    k = fitted_params(2);
    delta_1 = fitted_params(3);
    delta_2 = fitted_params(4);
    p = fitted_params(5);
    c = fitted_params(6);

    % Within-host ODE
    init_conditions = [167; 85; 69643];
    dydt = @(t,y)[lambda_w - k*y(1)*y(3) - delta_1*y(1);
                  k*y(1)*y(3) - delta_2*y(2);
                  p*y(2) - c*y(3)];
    [~, ysol] = ode23s(dydt, treatment_timeforward, init_conditions);
    Vvec = ysol(:,3);
    
   
    N_steps = length(tforward);
    M = length(treatment_timeforward);
    S = zeros(N_steps, 1); 
    I = zeros(N_steps, 1); 
    A = zeros(N_steps, 1); 
    D = zeros(N_steps, 1); 
    d_old = zeros(M,1); 
    d_new = zeros(M,1);
    S(1) = 820; 
    I(1) = 42; 
    A(1) = 0.337;
    D(1) = trapz(d_old)*d_tau; 
    Bd_Integral = zeros(N_steps, 1); 
    N_total = zeros(N_steps, 1); 
    N_total(1) = S(1) + I(1) + D(1);

    % Time-dependent \tilde{\beta}_d(t)
    beta_d_fun = @(t) ...
        (t <= t_start_transition).*beta_d0 + ...
        (t > t_start_transition & t <= t_end_transition).* ...
        ((beta_d_post - beta_d0)/(t_end_transition - t_start_transition).*(t - t_start_transition) + beta_d0) + ...
        (t > t_end_transition).*beta_d_post;

    % Simulation loop
    for kstep = 1:N_steps-1
        t_curr = tforward(kstep) + base_year;
        beta_d_tilde_now = beta_d_fun(t_curr);
        tind = min(kstep, M);
        beta_d_vec = beta_d_tilde_now * Vvec;
        
        N_total(kstep) = S(kstep) + I(kstep) + D(kstep);
        if N_total(kstep) == 0
            N_total(kstep) = 1;
        end
        
        Bd_Integral(kstep) = d_tau * trapz(beta_d_vec .* d_old) / N_total(kstep);

        % Update equations
        S(kstep+1) = (S(kstep) + dt*Lambda) / ...
            (1 + dt*(beta*I(kstep)/N_total(kstep) + 0.02));
        
        I(kstep+1) = (I(kstep) + dt*(beta*S(kstep+1)*I(kstep)/N_total(kstep) + ...
            S(kstep+1)*Bd_Integral(kstep))) / ...
            (1 + (0.02 + 0.002 + rho)*dt);
        
        d_new(1) = (1/kappa)*rho*I(kstep+1);
        d_new(2:end) = d_old(1:end-1) ./ (1 + dt*0.02);

        A(kstep+1) = (A(kstep) + dt*(gamma*I(kstep+1))) / (1 + dt*(0.02 + 0.002));
        D(kstep+1) = trapz(d_new)*d_tau;
        d_old = d_new;
    end

    HIV_Inc = beta*S(1:end-1).*I(1:end-1)./N_total(1:end-1) + S(1:end-1).*Bd_Integral(1:end-1);
    HIV_Diag = rho*I(1:end-1);
    years = tforward(1:end-1) + base_year;
end

function [years, R0_series, R0_base, R0_coef, R0_pre2019] = ...
    compute_R0_series(params, dt, T_end, base_year, beta_d0, beta_d_post, t_start_transition, t_end_transition)
  
    beta  = params(2);
    gamma = params(3);
    rho   = params(4);

    
    mu     = 0.02;     
    kappa  = 52;      
    d_tau  = 0.52;     
    Treat_Age = 20*52; 

    % Within-host ODE to get V(τ)  
    t_tau = (0:d_tau:Treat_Age).';  
    lambda_w = 33.84; k = 0.0087; delta_1 = 0.00036; delta_2 = 0.085; p = 0.45; c = 0.96;
    initWH = [167; 85; 69643];
    dydt = @(t,y)[lambda_w - k*y(1)*y(3) - delta_1*y(1); ...
                  k*y(1)*y(3) - delta_2*y(2); ...
                  p*y(2) - c*y(3)];
    [~, ysol] = ode23s(dydt, t_tau, initWH);
    Vtau = ysol(:,3);                    % Mx1

    
    C = trapz(t_tau, Vtau .* exp(-(mu/kappa).*t_tau));   
   
    denom       = (mu + gamma + rho);
    R0_base     = beta / denom;
    R0_coef     = (rho / (kappa * denom)) * C;          
    R0_pre2019  = R0_base + R0_coef * beta_d0;

    % Time axis (columns)
    tforward = (0:dt:T_end).';
    years    = tforward + base_year;

%     % Build \tilde{beta}_d(t) and R0(t) with a simple loop 
%     R0_series = nan(size(years));                          % column
%     for i = 1:numel(years)
%         t = years(i);
%         if t <= t_start_transition
%             beta_d_t = beta_d0;
%         elseif t <= t_end_transition
%             beta_d_t = (beta_d_post - beta_d0)/(t_end_transition - t_start_transition) ...
%                        * (t - t_start_transition) + beta_d0;
%         else
%             beta_d_t = beta_d_post;
%         end
%         R0_series(i) = R0_base + R0_coef * beta_d_t;      
%     end
% end



% Build \tilde{beta}_d(t) and R0(t) with two pieces: t <= 2019 vs t > 2019
R0_series = nan(size(years));   % column

slope = -(beta_d0) / (t_end_transition - t_start_transition);  
% Here: slope = -0.0012 / (2030 - 2019)

for i = 1:numel(years)
    t = years(i);

    if t <= t_start_transition
        % Before or at 2019: constant
        beta_d_t = beta_d0;
    else
        % After 2019: linear decline
        beta_d_t = beta_d0 + slope * (t - t_start_transition);
    end

    R0_series(i) = R0_base + R0_coef * beta_d_t;
end
end

function [c, ceq] = beta_constraints(beta_d_post, beta)
   
    c = beta_d_post - 0.1 * beta;  
    ceq = [];                     
end