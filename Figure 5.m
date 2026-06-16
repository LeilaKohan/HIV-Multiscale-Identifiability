close all; 
close all;
clc;

year_of_infection = [0, 1, 2, 3, 4, 5, 9, 15, 16, 18, 19, 25, 26, 27, 28, ...   
  29, 30, 31, 32, 33, 34, 35, 36, 37, 38]';

Estimated_HIV_incident_numbers = [20.77, 65.39, 65.13, 131.07, 131.37,...
    85.37, 85.15, 49.53, 49.89, 58.6, 59.19, 56.40, 56.49, 53.70,...
    42.21, 41.10, 40.28, 38.89, 38.07, 38.66, 38.95,...
    38.68, 37.85, 37.30, 36.19]'; % Number of HIV incidents in 1000s

year_of_infection_post = [39, 40, 41]';   % 2020–2022

Estimated_HIV_incident_numbers_post = [34.2, 32.7, 31.8]'; % in 1000s

years_diagnose = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]';

HIV_diagnose_cases = [46703, 44324, 42775, 40945, 40183, 38931, 39653,...
    39543, 39163, 38048, 37132, 36350]'/1000; % CDC DATA

years_diagnose_post = [39, 40, 41, 42, 43]'; % 2020–2024

HIV_diagnose_cases_post = [30.403, 35.763, 37.660, 38.793, 34.788]'; % in 1000s


years_AIDS_classification = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, ...
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,...
    34, 35, 36, 37, 38]';

AIDS_classidication_cases = [11808, 18806, 27929, 34432, 41588, 47387, 57222,...
    69709, 73971, 67469, 64353, 56500, 45946, 39343, 37728, 38049, 36622,...
    36429, 36992, 35900, 33956, 32616, 31821, 31234, 30063, 27328, 25556,...
    24608, 23574, 19244, 18539, 18284, 17727, 17106, 16410]'/1000; % CDC 

years_AIDS_death = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, ...
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,...
    34, 35, 36, 37, 38]';

AIDS_death_cases = [6712, 11727, 15703, 19980, 26438, 30051, 35145, 39861,...
    44158, 49396, 50468, 37450, 22056, 18857, 18587, 18275, 18397, 18339,...
    18202, 17799, 17470, 17094, 16288, 15361, 14847, 13703, 13330, 13041,...
    12977, 13060, 12779, 12945, 12825, 12649, 12471]'/1000;



d_tau = 0.52;                 
kappa = 52;                  
dt = d_tau / kappa;          

Time_Epidemic = 50;           % 1981-2030 
Treatmeant_Age = 20 * 52;     % within-host horizon (20 yr)

base_year = 1981;


lb = [0, 0, 0.02, 1/5, 0];
ub = [6000, 1000, 1, 2, 100];
params0 = [8.38663746616526, 0.306293687587148, 0.120967014736454,...
    0.200000005673198, 0.00119503317441228];

nonlcon = @constraints_betas;
options = optimset('Display','iter');
k = fmincon(@(k) err_in_data(k, dt), params0, [], [], [], [], lb, ub, nonlcon, options);

fprintf('  Lambda  = %g\n', k(1));
fprintf('  beta    = %g\n', k(2));
fprintf('  gamma   = %g\n', k(3));
fprintf('  rho     = %g\n', k(4));
fprintf('  beta_d  = %g\n', k(5));

% R0 (piecewise: constant ≤2019; changes after 2019) 
[R0_base, R0_diag, R0_pre2019, R0_post2019] = compute_R0_piecewise(k);

fprintf('  R0 (base term)          = %.6f\n', R0_base);
fprintf('  R0 (diagnosed term)     = %.6f\n', R0_diag);
fprintf('  R0 (≤ 2019, total)      = %.6f\n', R0_pre2019);
fprintf('  R0 (≥ 2020, total)      = %.6f\n', R0_post2019);


tforward = 0:dt:Time_Epidemic; 

[HIV_Cases, HIV_Diagnose, AIDS_Cases, AIDS_deaths, I, D] = ...
    FiniteDifference_HIV_Between_Host_Model(k, dt, Time_Epidemic);



real_t_inc   = tforward(1:end-1) + base_year;  
real_t_full  = tforward + base_year;            

% CDC-year helpers (add base_year)
real_year_inf     = year_of_infection + base_year;
real_year_diag    = years_diagnose    + base_year;
real_year_AIDScls = years_AIDS_classification + base_year;

%  time series of R0: constant up to 2019; post-2019 removes diagnosed transmission
R0_series = R0_post2019 * ones(size(real_t_full));
R0_series(real_t_full < 2019) = R0_pre2019;


%R0 figure 
figure('Position',[100,100,900,380]); hold on;

% color palette
cR0  = [0 0.4470 0.7410];  % blue
%cR0  = [0.8500 0.3250 0.0980];  % orange for R0(t)
cCP  = [0 0 0];                  % black for change point
cThr = [0 0.6 0];                % green for threshold

% curves + guide lines

hR0 = plot(real_t_full(:), R0_series(:), '-b', 'LineWidth', 2.2, 'Color', cR0);
hCP  = xline(2019, '--', 'LineWidth', 1.4, 'Color', cCP);
hThr = yline(1,    '--', 'LineWidth', 1.4, 'Color', cThr);


xlabel('Year', 'FontWeight','bold', 'FontSize',12);
ylabel('Basic reproduction number  R_0', 'FontWeight','bold', 'FontSize',12);


% title({'Refitting $\tilde{\beta}_d$ after 2019', ...
%        'Basic reproduction number $R_0(t)$'}, ...
%       'Interpreter','latex');
title({'Refitting $\tilde{\beta}_d$ after 2019', ...
       'Basic reproduction number $R_0(t)$'}, ...
      'Interpreter','latex', ...
      'FontWeight','bold');

lgd = legend([hR0 hCP hThr], ...
    {'$R_0(t)$', 'Change point: 2019', 'Threshold $R_0=1$'}, ...
    'Location','northeast', 'Interpreter','latex', 'Box','off');
set(lgd,'FontSize',13);

grid on; box on;
xlim([1981, 2032]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
hold off;




%  INCIDENCE FIGURE 1981-2030

figure('Position', [100, 100, 1000, 500]);

% 2019 vertical dashed line
xline(2019, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Goal triangle markers
scatter(2025, 9.3, 100, [0 0.6 0], 'v', 'filled', 'HandleVisibility', 'off');
scatter(2030, 3.0, 100, [0 0 0.8], 'v', 'filled', 'HandleVisibility', 'off');


% model line and CDC dots 
p1 = plot(real_t_inc, HIV_Cases, '-b', 'LineWidth', 2); hold on
p2 = scatter(real_year_inf, Estimated_HIV_incident_numbers, ...
             80, [0.8 0 0.8], 'filled');
% CDC ≥2020 
real_year_inf_post = year_of_infection_post + base_year;
p3 = scatter(real_year_inf_post, Estimated_HIV_incident_numbers_post, ...
             80, [0.85 0 0], 'filled');

% Legend
legend([p1 p2 p3], ...
   {'Model predicted HIV incidences', ...
    'CDC reported incidences (≤2019)', ...
    'CDC reported incidences (2020+)'}, ...
   'Location', 'northwest', 'AutoUpdate','off');
 
xline(2019,'--k','LineWidth',1.2,'HandleVisibility','off');
scatter(2025,9.3,100,[0 0.6 0],'v','filled','HandleVisibility','off');
scatter(2030,3.0,100,[0 0 0.8],'v','filled','HandleVisibility','off');

xlabel('Years of HIV Epidemic in USA', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('HIV Cases Per Year (in 1000s)', 'FontWeight', 'bold', 'FontSize', 12);


title('HIV Incidences Projected to 2030');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;

xlim([1981, 2032]);

% Dashed vertical line at 2019 
yl = ylim;
xline(2019, '--k', 'LineWidth', 1.2);
text(2013 + 0.2, yl(2) * 0.75, ...
    'Ending the HIV Epidemic Initiated', ...
    'FontSize', 14, 'Color', [0 0.5 0]);

% Plot 2025 and 2030 goal markers (triangles)
goal2025 = 9.3;   % in 1000s
goal2030 = 3.0;

scatter(2025, goal2025, 100, [0 0.6 0], 'v', 'filled');  
scatter(2030, goal2030, 100, [0 0 0.8], 'v', 'filled');  

% arrow labels 
ax = gca;
ax_pos = get(ax, 'Position');         
xlim_ax = xlim(ax);
ylim_ax = ylim(ax);

% Converters from data to normalized figure units
x2fig = @(x) ax_pos(1) + (x - xlim_ax(1)) / diff(xlim_ax) * ax_pos(3);
y2fig = @(y) ax_pos(2) + (y - ylim_ax(1)) / diff(ylim_ax) * ax_pos(4);

arrow_offset = 40;  % vertical arrows (in 1000s)

% 2025 annotation (green)
annotation('textarrow', ...
    [x2fig(2025 - 2.5), x2fig(2025)], ...
    [y2fig(goal2025 + arrow_offset), y2fig(goal2025)], ...
    'String', '2025 Year Target: 9.3k', ...
    'FontSize', 12, 'Color', [0 0.5 0]);

% 2030 annotation (blue)
annotation('textarrow', ...
    [x2fig(2030 - 3), x2fig(2030)], ...
    [y2fig(goal2030 + arrow_offset), y2fig(goal2030)], ...
    'String', '2030 Year Target: 3k', ...
    'FontSize', 12, 'Color', [0 0 0.6]);

hold off;


goal2025_diag = 9.588;  % in 1000s
goal2030_diag = 3.0;

figure('Position', [100, 100, 1000, 500]);

% Plot model predictions
p1 = plot(real_t_full, HIV_Diagnose, '-b', 'LineWidth', 2); hold on;

% CDC observed data
p2a = scatter(real_year_diag, HIV_diagnose_cases, ...
             80, [0.8 0 0.8], 'filled');

p3 = scatter(real_year_diag_post, HIV_diagnose_cases_post, ...
              80, 'o', 'MarkerFaceColor',[0.85 0 0], 'MarkerEdgeColor','k');

legend([p1 p2a p3], ...
   {'Model predicted HIV diagnoses', ...
    'CDC reported diagnoses (≤2019)', ...
    'CDC reported diagnoses (2020+)'}, ...
   'Location','northwest','AutoUpdate','off');


% 2019 vertical dashed line
xline(2019, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Goal triangle markers
scatter(2025, goal2025_diag, 100, [0 0.6 0], 'v', 'filled', 'HandleVisibility', 'off');
scatter(2030, goal2030_diag, 100, [0 0 0.8], 'v', 'filled', 'HandleVisibility', 'off');

% 
% legend([p1 p2], {'Model Predicted HIV Diagnoses', 'CDC Reported Diagnoses'}, ...
%        'Location', 'northwest', ...
%        'AutoUpdate','off');

xlabel('Years of HIV Epidemic in USA', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('HIV Diagnoses Per Year (in 1000s)', 'FontWeight', 'bold', 'FontSize', 12);
title('HIV Diagnoses Projected to 2030');
grid on;
xlim([1981, 2032]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

yl = ylim;
text(2013 + 0.2, yl(2) * 0.75, ...
     'Ending the HIV Epidemic Initiated', ...
     'FontSize', 14, 'Color', [0 0.5 0]);

% Add triangles
scatter(2025, goal2025_diag, 100, [0 0.6 0], 'v', 'filled');
scatter(2030, goal2030_diag, 100, [0 0 0.8], 'v', 'filled');

% arrow annotations
ax = gca;
ax_pos = get(ax, 'Position');
xlim_ax = xlim(ax); ylim_ax = ylim(ax);
x2fig = @(x) ax_pos(1) + (x - xlim_ax(1)) / diff(xlim_ax) * ax_pos(3);
y2fig = @(y) ax_pos(2) + (y - ylim_ax(1)) / diff(ylim_ax) * ax_pos(4);
arrow_offset = 40;

% 2025 Arrow
annotation('textarrow', ...
    [x2fig(2025 - 2.5), x2fig(2025)], ...
    [y2fig(goal2025_diag + arrow_offset), y2fig(goal2025_diag)], ...
    'String', '2025 Year Target: 9.6k', ...
    'FontSize', 12, 'Color', [0 0.5 0]);

% 2030 Arrow
annotation('textarrow', ...
    [x2fig(2030 - 3), x2fig(2030)], ...
    [y2fig(goal2030_diag + arrow_offset), y2fig(goal2030_diag)], ...
    'String', '2030 Year Target: 3k', ...
    'FontSize', 12, 'Color', [0 0 0.6]);

hold off;


% AIDS-classification panel
figure('Position', [100, 100, 1000, 500]);

plot(real_t_inc, AIDS_Cases, '-b','LineWidth',2); hold on
scatter(real_year_AIDScls, AIDS_classidication_cases, ...
        80, [0.8 0 0.8], 'filled')

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('AIDS Classifications (in 1000s)', 'FontWeight', 'bold', 'FontSize', 12);
title('Annual AIDS Classifications in the U.S.');
grid on;
set(gca,'FontSize',12,'LineWidth',1.4);

xlim([1981, 2032]);  
hold off;


% Data prep
ID_sum = I + D;
ID_sum(ID_sum == 0) = eps;
D_over_ID = D ./ ID_sum;

real_t = (0:length(D)-1) * dt + base_year;

cdc_years     = [2017, 2018, 2019, 2020, 2021, 2022];
cdc_knowledge = [85.8, 85.9, 86.3, 86.5, 86.8, 87.2] / 100;


figure('Color','w','Position',[100 100 1000 500]); hold on;

% Model-predicted line
plot(real_t, D_over_ID, '-b', 'LineWidth', 2);

% CDC observed data (filled circles)
scatter(cdc_years, cdc_knowledge, 80, [0.8 0 0.8], 'filled');


% Draw the horizontal line without label
yline(0.95, '--', 'Color', [0 0.6 0], 'LineWidth', 1.2);

% Manually place the label text
text(1987, 0.92, '2025 & 2030 Goal: 95%', ...
     'FontSize', 12, ...
     'Color', [0 0.6 0], ...
     'FontWeight', 'normal');


% Dashed vertical line at 2019 
yl = ylim;
xline(2019, '--k', 'LineWidth', 1.2);
text(2013 + 0.2, yl(2) * 0.6, ...
    'Ending the HIV Epidemic Initiated', ...
    'FontSize', 14, 'Color', [0 0.5 0]);

xlabel('Years of HIV Epidemic in USA', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Proportion with Knowledge of HIV Status', 'FontWeight', 'bold', 'FontSize', 12);
title('Knowledge of HIV Status: Model vs. CDC & National Goal');

legend({'Model predicted D / (I + D)', ...
        'CDC reported knowledge'}, ...
        'Location','southeast', 'FontSize', 11);

xlim([1981 2032]);


set(gca, 'FontSize', 12, ...
         'LineWidth', 1.5, ...
         'XGrid','on', 'YGrid','on', ...
         'Box','on');

hold off;


%Percentage of people with HIV who have viral suppression
% years = 2017:2022;
% suppression = [0.631, 0.647, 0.655, 0.646, 0.659, 0.651];  
% 
% figure;
% plot(years, suppression, 'o-', 'LineWidth',2);
% %yline(0.95, '--g', '95% Goal');
% xlabel('Year'); ylabel('Viral Suppression (%)');
% title('Viral Suppression Among Diagnosed PLWH');
% ylim([0.5 1]); 
% grid on;



function [HIV_Cases, HIV_Diagnose, AIDS_Cases, AIDS_deaths, I, D] ...
          = FiniteDifference_HIV_Between_Host_Model(params, dt, Time_Epidemic)

Lambda = params(1); 
beta  = params(2); 
gamma = params(3);
rho    = params(4); 
beta_d_coeff = params(5);
mu = 0.02; 
alpha = 0.002;


d_tau = 0.52; 
kappa = 52;  
Treatmeant_Age = 20*52;
T0=167; 
T_i0=85; 
V0=69643;

[t_sol, y_sol] = ode23s(@(t,y) HIV_Within_Host(y), ...
                        0:d_tau:Treatmeant_Age, [T0; T_i0; V0]);
     
% After generating beta_d_vec
% t_within = 0:d_tau:Treatmeant_Age;
base_year = 1981;
beta_d_vec = beta_d_coeff * y_sol(:,3);

% Force beta_d to zero after 2019 for prediction only
% year_vec = base_year + t_within;
%beta_d_vec(year_vec >= 2019 + dt) = 0;


tforward = 0:dt:Time_Epidemic;  
N = numel(tforward);
S = zeros(N,1); 
I = zeros(N,1); 
A = zeros(N,1);
d_old = zeros(numel(beta_d_vec),1); 
d_new = d_old;
D = zeros(N,1); 
Bd_Int = zeros(N,1);  Ntot = zeros(N,1);

S(1)=820; 
I(1)=42; 
A(1)=0.337;
D(1)=0; 
Ntot(1)=S(1)+I(1)+D(1);
 


    fade_cutoff = 2019;

for k = 1:N-1
    Ntot(k) = S(k) + I(k) + D(k);
    year_k = base_year + (k - 1) * dt;

    if Ntot(k) == 0
        Ntot(k) = eps;  % Avoid division by zero
    end

    % Compute Bd_Int only before fade_cutoff
    if year_k <= fade_cutoff
        Bd_Int(k) = d_tau * trapz(beta_d_vec .* d_old) / Ntot(k);
    else
        Bd_Int(k) = 0;
    end

    % Update state variables using finite difference
    S(k+1) = (S(k) + dt * Lambda) / (1 + dt * (beta * I(k) / Ntot(k) + mu));

    I(k+1) = (I(k) + dt * (beta * S(k+1) * I(k) / Ntot(k) + ...
              S(k+1) * Bd_Int(k))) / ...
             (1 + dt * (mu + alpha + rho));

    d_new(1) = (1/kappa) * rho * I(k+1);
    d_new(2:end) = d_old(1:end-1) / (1 + dt * mu);

    A(k+1) = (A(k) + dt * gamma * I(k+1)) / (1 + dt * (mu + alpha));
    D(k+1) = trapz(d_new) * d_tau;

    d_old = d_new;
end


Ntot(end) = S(end) + I(end) + D(end);

HIV_Cases    = beta * S(1:end-1) .* I(1:end-1) ./ Ntot(1:end-1) + ...
               S(1:end-1) .* Bd_Int(1:end-1);
HIV_Diagnose = rho * I;
AIDS_Cases   = gamma * I(1:end-1);
AIDS_deaths  = (mu + alpha) * A;

end


function dy = HIV_Within_Host(y)

lambda=33.84;
k=0.0087; 
delta1=0.00036; 
delta2=0.085;
p=0.45; 
c=0.96;
T=y(1); 
Ti=y(2); 
V=y(3);
dy=zeros(3,1);
dy(1) = lambda - k*T*V - delta1*T;
dy(2) = k*T*V - delta2*Ti;
dy(3) = p*Ti - c*V;
end

function err = err_in_data(k,dt)
[HIV_Inc, HIV_Diag, AIDS_Cases, ~] = FiniteDifference_HIV_Between_Host_Model(k, dt, 49);

year_inc_idx  = round([0,1,2,3,4,5,9,15,16,18,19,25,26,27,28,29,30,31,32,...
                       33,34,35,36,37,38]/dt) + 1;
year_diag_idx = round((27:38)/dt) + 1;
year_cls_idx  = round((4:38)/dt)  + 1;

CDC_inc  = [20.77,65.39,65.13,131.07,131.37,85.37,85.15,49.53,49.89, ...
            58.6,59.19,56.40,56.49,53.70,42.21,41.10,40.28,38.89, ...
            38.07,38.66,38.95,38.68,37.85,37.30,36.19]';
CDC_diag = [46703,44324,42775,40945,40183,38931,39653,39543,39163, ...
            38048,37132,36350]' / 1000;
CDC_cls  = [11808,18806,27929,34432,41588,47387,57222,69709,73971, ...
            67469,64353,56500,45946,39343,37728,38049,36622,36429, ...
            36992,35900,33956,32616,31821,31234,30063,27328,25556, ...
            24608,23574,19244,18539,18284,17727,17106,16410]' / 1000;

err = 2*sum((HIV_Inc(year_inc_idx)  - CDC_inc ).^2) + ...
      20*sum((HIV_Diag(year_diag_idx) - CDC_diag).^2) + ...
      11*sum((AIDS_Cases(year_cls_idx) - CDC_cls).^2);
end


function [R0_base, R0_diag, R0_pre, R0_post] = compute_R0_piecewise(params)

   
    beta         = params(2);
    gamma        = params(3);
    rho          = params(4);
    beta_d_coeff = params(5);

 
    mu    = 0.02;      
    kappa = 52;       
    d_tau = 0.52;     
    Treatmeant_Age = 20*52;  

    % Within-host
    T0 = 167; Ti0 = 85; V0 = 69643;

    % Get V(τ) to build beta_d(τ)
    [tau, y] = ode23s(@(t,y) HIV_Within_Host(y), 0:d_tau:Treatmeant_Age, [T0; Ti0; V0]);
    V_tau      = y(:,3);
    beta_d_tau = beta_d_coeff .* V_tau;

    % gamma_d(τ) = 0  ⇒  ∫_0^τ (mu + gamma_d) ds = mu * τ
    exp_term = exp(-(mu/kappa) .* tau);

    % Outer integral
    integral_term = trapz(tau, beta_d_tau .* exp_term);

    % Components
    denom   = (mu + gamma + rho);
    R0_base = beta / denom;
    R0_diag = (rho / (kappa * denom)) * integral_term;

    % Piecewise totals
    R0_pre  = R0_base + R0_diag;  % ≤ 2019 (includes diagnosed transmission)
    R0_post = R0_base;            % ≥ 2020 (diagnosed transmission removed)
end


function [c, ceq] = constraints_betas(p)
c   = p(5) - 0.1*p(2); 
ceq = [];
end

