clear all
close all
clc

num_iter = 5000;

year_of_infection = [0, 1, 2, 3, 4, 5, 9, 15, 16, 18, 19, 25, 26, 27, 28, ...   
  29, 30, 31, 32, 33, 34, 35, 36, 37, 38]';

Estimated_HIV_incident_numbers = [20.77, 65.39, 65.13, 131.07, 131.37,...
    85.37, 85.15, 49.53, 49.89, 58.6, 59.19, 56.40, 56.49, 53.70,...
    42.21, 41.10, 40.28, 38.89, 38.07, 38.66, 38.95,...
    38.68, 37.85, 37.30, 36.19]'; % Number of HIV incidents in 1000s

years_diagnose = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]';

HIV_diagnose_cases = [46703, 44324, 42775, 40945, 40183, 38931, 39653,...
    39543, 39163, 38048, 37132, 36350]'/1000; % CDC DATA


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
    kappa = 52;  % Scaling factor for different time scales

    dt = d_tau/kappa; 

    Time_Epidemic = 39;   %1981-2025 44 years
    Treatmeant_Age = 20*52;   %max time in 80 patients 4years but treatment take for almost 20 years

    tforward = 0:dt:Time_Epidemic;


   %params = [Lambda  beta gamma rho betad_tilde];

 lb = [0, 0, 0.02, 1/5, 0];
 ub = [6000, 1000, 1, 2, 100];


   params = [8.38663746616526	0.306293687587148	0.120967014736454	0.200000005673198	0.00119503317441228];
 


year_infection_measure = [0, 1, 2, 3, 4, 5, 9, 15, 16, 18, 19, 25, 26, 27, 28, ...   
  29, 30, 31, 32, 33, 34, 35, 36, 37, 38]/dt+1;

years_diagnose_measure = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]/dt+1;

years_AIDS_classification_measure = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, ...
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,...
    34, 35, 36, 37, 38]/dt+1;

[HIV_Cases, HIV_Diagnose, AIDS_Cases, AIDS_deaths, I, D] = FiniteDifference_HIV_Between_Host_Model(params);


Model_HIV_Inc = HIV_Cases(year_infection_measure);
Model_HIV_Diagnoses = HIV_Diagnose(years_diagnose_measure);
Model_AIDS_Cases = AIDS_Cases(years_AIDS_classification_measure);


res_inc  = Model_HIV_Inc  - Estimated_HIV_incident_numbers;
res_diag = Model_HIV_Diagnoses - HIV_diagnose_cases;
res_aids = Model_AIDS_Cases - AIDS_classidication_cases;
sd_inc  = std(res_inc ./ Model_HIV_Inc);
sd_diag = std(res_diag ./ Model_HIV_Diagnoses);
sd_aids = std(res_aids ./ Model_AIDS_Cases);

HIV_Inc_Data_Set = zeros(length(Model_HIV_Inc),num_iter);
HIV_Diagnoses_Data_Set = zeros(length(Model_HIV_Diagnoses),num_iter);
AIDS_Data_Set = zeros(length(Model_AIDS_Cases),num_iter);
estimated_params = zeros(num_iter, length(params));


    parfor iter = 1:num_iter
       
        HIV_Inc_noisy = Model_HIV_Inc + sd_inc*randn(size(Model_HIV_Inc)).*Model_HIV_Inc;
        HIV_Diagnoses_noisy = Model_HIV_Diagnoses + sd_diag*randn(size(Model_HIV_Diagnoses)).*Model_HIV_Diagnoses;
        AIDS_Cases_noisy = Model_AIDS_Cases + sd_aids*randn(size(Model_AIDS_Cases)).*Model_AIDS_Cases;
 

        HIV_Inc_Data_Set(:,iter) = HIV_Inc_noisy;
        HIV_Diagnoses_Data_Set(:,iter) = HIV_Diagnoses_noisy;
        AIDS_Data_Set(:,iter) = AIDS_Cases_noisy;
        
        err_fun = @(k)err_in_data_with_noisy(k, HIV_Inc_noisy, HIV_Diagnoses_noisy, AIDS_Cases_noisy);


        
        k_est = fmincon(err_fun, params, [], [], [], [], lb, ub, @constraints_betas);
        estimated_params(iter,:) = k_est;

    end

HIV_Inc_Model_Pred = zeros(length(HIV_Cases),num_iter);
HIV_Diagnoses_Model_Pred = zeros(length(HIV_Diagnose),num_iter);
AIDS_Model_Pred = zeros(length(AIDS_Cases),num_iter);

 for iter=1:num_iter
     
     k_params = estimated_params(iter,:);
     [HIV_Cases, HIV_Diagnose, AIDS_Cases, ~, ~, ~] = FiniteDifference_HIV_Between_Host_Model(k_params );
     HIV_Inc_Model_Pred(:,iter) = HIV_Cases;
     HIV_Diagnoses_Model_Pred(:,iter) = HIV_Diagnose;
     AIDS_Model_Pred(:,iter) = AIDS_Cases;
 end

base_year = 1981;
real_tforward = tforward + base_year;
real_tforward_trimmed = tforward(1:end-1) + base_year;
real_year_of_infection = year_of_infection + base_year;
real_years_diagnose = years_diagnose + base_year;
real_years_AIDS_classification = years_AIDS_classification + base_year;
real_years_AIDS_death = years_AIDS_death + base_year;


% Photos

%HIV Incidence Over Time
figure;
plot(real_tforward_trimmed, HIV_Cases, '-b', 'LineWidth', 2); hold on;
plot(real_year_of_infection, Estimated_HIV_incident_numbers, 'r.', 'MarkerSize', 25);
for i=1:num_iter
    plot(real_year_of_infection, HIV_Inc_Data_Set(:,i), '.', 'MarkerSize', 25);
    plot(real_tforward_trimmed, HIV_Inc_Model_Pred(:,i), '-b', 'LineWidth', 2.5);
end
xlabel('Year');
ylabel('HIV Incidence (in thousands)');
title('(a) Estimated Annual HIV Incidence in the U.S.');
%legend('Model', 'CDC Data', 'Location', 'NorthEast');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);


% Diagnosed HIV Population Over Time
figure;
plot(real_tforward, HIV_Diagnose, '-b', 'LineWidth', 2); hold on;
plot(real_years_diagnose, HIV_diagnose_cases, 'r.', 'MarkerSize', 25);
for i=1:num_iter
    plot(real_years_diagnose, HIV_Diagnoses_Data_Set(:,i), '.', 'MarkerSize', 25);
    plot(real_tforward, HIV_Diagnoses_Model_Pred(:,i), '-b', 'LineWidth', 2)
end
xlabel('Year');
ylabel('HIV Diagnoses (in thousands)');
title('(b) Annual HIV Diagnoses in the U.S.');
%legend('Model', 'CDC Data', 'Location', 'NorthEast');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5); 


%AIDS Classifications Over Time
figure;
plot(real_tforward_trimmed, AIDS_Cases, '-b', 'LineWidth', 2); hold on;
plot(real_years_AIDS_classification, AIDS_classidication_cases, 'r.', 'MarkerSize', 25);

for i=1:num_iter
    plot(real_years_AIDS_classification, AIDS_Data_Set(:,i), '.', 'MarkerSize', 25);
    plot(real_tforward_trimmed, AIDS_Model_Pred(:,i), '-b', 'LineWidth', 0.2);
end
xlabel('Year');
ylabel('AIDS Classifications (in thousands)');
title('(c) Annual AIDS Classifications in the U.S.');
%legend('Model', 'CDC Data', 'Location', 'NorthEast');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);




%%%%%%%

z99 = 2.576;   % 99% two-sided
z3  = 3.0;     % μ ± 3σ (~99.7%)
lower_bound_Inc = HIV_Cases.* (1 - 2.576*sd_inc);
upper_bound_Inc = HIV_Cases.* (1 + 2.576*sd_inc);

lower_bound_diag = HIV_Diagnose.* (1 - 2.576*sd_diag);
upper_bound_diag = HIV_Diagnose.* (1 + 2.576*sd_diag);

lower_bound_aids = AIDS_Cases.* (1 - 2.576*sd_aids);
upper_bound_aids = AIDS_Cases.* (1 + 2.576*sd_aids);



% === HIV Incidence Over Time ===
Estimated_HIV_incident_numbers(Estimated_HIV_incident_numbers < 0) = NaN;
HIV_Inc_Data_Set(HIV_Inc_Data_Set < 0) = NaN;
lower_bound_Inc = max(lower_bound_Inc, 0);
figure('Color','w'); hold on;

% Virtual epidemic fits 
for i = 1:num_iter
    plot(real_tforward_trimmed, HIV_Inc_Model_Pred(:,i), '-', ...
        'Color', [0 0.6 0.7 0.08], 'LineWidth', 0.3, 'HandleVisibility','off');
end

% Simulated HIV incidence (magenta points)
for i = 1:num_iter
    plot(real_year_of_infection, HIV_Inc_Data_Set(:,i), '.', ...
        'MarkerSize', 5, 'Color', [0.8 0.2 0.8 0.1], 'HandleVisibility','off');
end

% Main model fit (thick navy line)
plot(real_tforward_trimmed, HIV_Cases, '-', ...
     'Color', [0 0.2 0.5], 'LineWidth', 3);

% CDC data 
plot(real_year_of_infection, Estimated_HIV_incident_numbers, 'o', ...
     'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 5);

xlabel('Year');
ylabel('HIV Incidence (in thousands)');
title('Estimated Annual HIV Incidence in the U.S.');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'GridAlpha', 0.3);


box on;
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.LineWidth = 1.2;

yl = ylim; xl = xlim;
plot(xl, [yl(2) yl(2)], 'k', 'LineWidth', 1.2); 
plot([xl(1) xl(1)], yl, 'k', 'LineWidth', 1.2); 

h1 = plot(NaN, NaN, '-', 'Color', [0 0.6 0.7], 'LineWidth', 1.2);       % teal
h2 = plot(NaN, NaN, '-', 'Color', [0.8 0.2 0.8], 'LineWidth', 1.2);     % magenta
h3 = plot(NaN, NaN, '-', 'Color', [0 0.2 0.5], 'LineWidth', 3);         % navy
h4 = plot(NaN, NaN, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0], 'MarkerSize', 6);

legend([h1,h2,h3,h4], ...
       {'Virtual Epidemic Fits', 'Simulated HIV Incidence', 'Model Fit to CDC', 'CDC Data'}, ...
       'Location', 'northeast', 'FontSize', 10, 'Box', 'off');




%  Annual HIV Diagnoses 
HIV_diagnose_cases(HIV_diagnose_cases < 0) = NaN;
HIV_Diagnoses_Data_Set(HIV_Diagnoses_Data_Set < 0) = NaN;

figure('Color','w'); hold on;

%  Virtual epidemic fits 
for i = 1:num_iter
    plot(real_tforward, HIV_Diagnoses_Model_Pred(:,i), '-', ...
        'Color',[0 0.6 0.7 0.08], 'LineWidth', 0.3, 'HandleVisibility','off');
end

%  Simulated diagnoses (magenta points)
for i = 1:num_iter
    plot(real_years_diagnose, HIV_Diagnoses_Data_Set(:,i), '.', ...
        'MarkerSize', 5, 'Color', [0.8 0.2 0.8 0.1], 'HandleVisibility','off');
end

%  Main model fit 
plot(real_tforward, HIV_Diagnose, '-', ...
     'Color', [0 0.2 0.5], 'LineWidth', 3);

%  CDC data (red dots, black edge) 
plot(real_years_diagnose, HIV_diagnose_cases, 'o', ...
      'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 5);

xlabel('Year');
ylabel('HIV Diagnoses (in thousands)');
title('Annual HIV Diagnoses in the U.S.');
grid on; set(gca,'FontSize',12,'LineWidth',1.2,'GridAlpha',0.3);


box on;
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.LineWidth = 1.2;
yl = ylim; xl = xlim;
plot(xl, [yl(2) yl(2)], 'k', 'LineWidth', 1.2); 
plot([xl(1) xl(1)], yl, 'k', 'LineWidth', 1.2); 


h1 = plot(NaN, NaN, '-', 'Color', [0 0.6 0.7], 'LineWidth', 1.2);       % teal
h2 = plot(NaN, NaN, '-', 'Color', [0.8 0.2 0.8], 'LineWidth', 1.2);     % magenta
h3 = plot(NaN, NaN, '-', 'Color', [0 0.2 0.5], 'LineWidth', 3);         % navy
h4 = plot(NaN, NaN, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0], 'MarkerSize', 6);
legend([h1,h2,h3,h4], ...
       {'Virtual Epidemic Fits', 'Simulated HIV Diagnoses', 'Model Fit to CDC', 'CDC Data'}, ...
       'Location', 'northeast', 'FontSize', 10, 'Box', 'off');


% Annual AIDS Classifications 
AIDS_classidication_cases(AIDS_classidication_cases < 0) = NaN;
AIDS_Data_Set(AIDS_Data_Set < 0) = NaN;

figure('Color','w'); hold on;

%  Virtual epidemic fits 
for i = 1:num_iter
    plot(real_tforward_trimmed, AIDS_Model_Pred(:,i), '-', ...
        'Color', [0 0.6 0.7 0.08], 'LineWidth', 0.3, 'HandleVisibility','off');
end

% Simulated AIDS classifications (magenta points)
for i = 1:num_iter
    plot(real_years_AIDS_classification, AIDS_Data_Set(:,i), '.', ...
        'MarkerSize', 5, 'Color', [0.8 0.2 0.8 0.10], 'HandleVisibility','off');
end

%  Main model fit 
plot(real_tforward_trimmed, AIDS_Cases, '-', ...
     'Color', [0 0.2 0.5], 'LineWidth', 3);

%  CDC data 
plot(real_years_AIDS_classification, AIDS_classidication_cases, 'o', ...
     'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0], 'MarkerSize', 5);

xlabel('Year');
ylabel('AIDS Classifications (in thousands)');
title('Annual AIDS Classifications in the U.S.');
grid on; set(gca,'FontSize',12,'LineWidth',1.2,'GridAlpha',0.3);


box on;
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.LineWidth = 1.2;
yl = ylim; xl = xlim;
plot(xl, [yl(2) yl(2)], 'k', 'LineWidth', 1.2); 
plot([xl(1) xl(1)], yl, 'k', 'LineWidth', 1.2); 


h1 = plot(NaN, NaN, '-', 'Color', [0 0.6 0.7], 'LineWidth', 1.2);       % teal
h2 = plot(NaN, NaN, '-', 'Color', [0.8 0.2 0.8], 'LineWidth', 1.2);     % magenta
h3 = plot(NaN, NaN, '-', 'Color', [0 0.2 0.5], 'LineWidth', 3);         % navy
h4 = plot(NaN, NaN, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0], 'MarkerSize', 6);
legend([h1,h2,h3,h4], ...
       {'Virtual Epidemic Fits','Simulated AIDS Classifications','Model Fit to CDC','CDC Data'}, ...
       'Location','northeast','FontSize',10,'Box','off');




function [HIV_Cases, HIV_Diagnose, AIDS_Cases, AIDS_deaths, I, D] = FiniteDifference_HIV_Between_Host_Model(params)

    d_tau = 0.52;
    kappa = 52;  % Scaling factor for different time scales

    dt = d_tau/kappa; 

    Time_Epidemic = 39;   %1981-2025 44 years
    Treatmeant_Age = 20*52;   %max time in 80 patients 4years but treatment take for almost 20 years

    % Parameters
    Lambda = params(1);       % Birth rate
    mu = 0.02;         % Natural death rate is fixed
    beta = params(2);       % Transmission rate
    gamma = params(3);      % AIDS progression rate
    rho = params(4);         % Treatment rate 
  
  
    alpha = 0.002;
    
    tforward = 0:dt:Time_Epidemic;
    treatment_timeforward = 0:d_tau:Treatmeant_Age;
    
    %treatment initial conditions
    T0 = 167;   
    T_i0 = 85;   
    V0 = 69643; 
    initial_condition = [T0; T_i0; V0];

    [t_sol, y_sol] = ode23s(@(t,y)HIV_Within_Host(y),treatment_timeforward,initial_condition);

  
    beta_d_vec = params(5)*y_sol(:,3); 

    N_steps = length(tforward);         % Number of time 
    M = length(treatment_timeforward);  % Number of treatment age 


    S = zeros(N_steps, 1);  
    I = zeros(N_steps, 1);  
    A = zeros(N_steps, 1);  

    d_old = zeros(M,1);
    d_new = zeros(M,1);

    D = zeros(N_steps, 1); 
 
    
    S(1) = 820;      
    I(1) = 42;    
    A(1) = 0.337;         
    d(1, :) = 0;      
 
    D(1) = trapz(d_old) * d_tau;
    
    Bd_Integral = zeros(N_steps, 1);
   
    N_total = zeros(N_steps, 1);

    N_total(1) = S(1) + I(1) + D(1);



    for k = 1:N_steps-1   
        
        N_total(k) = S(k) + I(k) + D(k);

       
        Bd_Integral(k) = d_tau*trapz(beta_d_vec.*d_old)/N_total(k);
        
        
        S(k+1) = (S(k) + dt * Lambda) / (1 + dt * (beta * I(k) / N_total(k) + ...
           mu));
        
        
        I(k+1) = (I(k) + dt * (beta * S(k+1) * I(k) / N_total(k) + ...
            S(k+1) * Bd_Integral(k))) / ...
            (1 + (mu + alpha + rho) * dt);
        

        d_new(1) = (1/kappa)*rho * I(k+1);
        d_new(2:end) = d_old(1:end-1)./(1 + dt*mu);
   
      
        
        A(k+1) = (A(k) + dt * (gamma * I(k+1))) / ...
            (1 + dt * (mu + alpha));

        D(k+1) = trapz(d_new) * d_tau; 

        d_old = d_new;

    end
    
HIV_Cases = beta*S(1:end-1).*I(1:end-1)./N_total(1:end-1) + S(1:end-1).*Bd_Integral(1:end-1);
HIV_Diagnose = rho*I;
AIDS_Cases = gamma*I(1:end-1);
AIDS_deaths = (mu + alpha)*A;

end

function error_in_data = err_in_data_with_noisy(k, HIV_Inc_noisy, HIV_Diagnoses_noisy, AIDS_Cases_noisy)


dt = 0.01;

year_infection_measure = [0, 1, 2, 3, 4, 5, 9, 15, 16, 18, 19, 25, 26, 27, 28, ...   
  29, 30, 31, 32, 33, 34, 35, 36, 37, 38]/dt+1;

years_diagnose_measure = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]/dt+1;

years_AIDS_classification_measure = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, ...
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,...
    34, 35, 36, 37, 38]/dt+1;



[HIV_Inc, HIV_Diagnose, AIDS_Cases, ~] = FiniteDifference_HIV_Between_Host_Model(k);

Model_HIV_Inc = HIV_Inc(year_infection_measure);
Model_HIV_Diagnoses = HIV_Diagnose(years_diagnose_measure);
Model_AIDS_Cases = AIDS_Cases(years_AIDS_classification_measure);



error_in_data = 2*sum((Model_HIV_Inc -  HIV_Inc_noisy).^2)+...
                20*sum((Model_HIV_Diagnoses - HIV_Diagnoses_noisy).^2) +...
                11*sum((Model_AIDS_Cases - AIDS_Cases_noisy).^2);
              


              
end

function [c,ceq] = constraints_betas(params)
ceq = [];
c = [params(5) - 0.1*params(2)];
end

    function dy = HIV_Within_Host(y) %inhost ode

      dy = zeros(3,1);
     


      fitted_params = [33.84, 0.0087, 0.00036, 0.085, 0.45, 0.96]; 

    lambda = fitted_params(1);   %target cell production rate
    k = fitted_params(2);         %infection rate
    delta_1 = fitted_params(3);  %death rate of target cells
    delta_2 = fitted_params(4);     %death rate of infected cells
    p = fitted_params(5);     %production rate of virus     
    c = fitted_params(6);   %clearence rate of virus

    T = y(1);
    T_i = y(2);
    V = y(3);


    dy(1) = lambda - k*T*V - delta_1*T;
    dy(2) = k*T*V - delta_2*T_i;
    dy(3) = p*T_i - c*V;

end


% === Bootstrap Histograms: 99% Trimmed Distributions ===
% Each column of estimated_params corresponds to: [Lambda, beta, gamma, rho, beta_tilde_d]


param_names = {'\Lambda', '\beta', '\gamma', '\rho', '\tilde{\beta}_d'};
num_params  = numel(param_names);

figure('Name','Histograms of Trimmed Bootstrap Distributions', ...
       'Position',[100 100 900 600]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for i = 1:num_params
    p = all_params(:,i);
    p = p(isfinite(p));  % clean NaNs/Infs

    % ---- Central 99% interval ----
    bounds = prctile(p,[0.5 99.5]);
    p_trim = p(p >= bounds(1) & p <= bounds(2));

    nexttile; hold on; box on; grid on


    % 3. Histogram of central 99% (blue bars)
    histogram(p_trim,50,'Normalization','probability',...
        'FaceColor',[0 0.447 0.741],'EdgeColor','k');


    xlabel(['$',param_names{i},'$'],'Interpreter','latex','FontSize',12);
    ylabel('Probability');
end


param_names = {'\Lambda', '\beta', '\gamma', '\rho', '\tilde{\beta}_d'};

% Virtual epidemic parameter estimates
true_vals = [8.39, 0.306, 0.121, 0.20, 0.0012];

num_params = numel(param_names);

figure('Name','Histograms of Trimmed Bootstrap Distributions', ...
       'Position',[100 100 900 600]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for i = 1:num_params
    p = all_params(:,i);
    p = p(isfinite(p));

    bounds = prctile(p,[0.5 99.5]);
    p_trim = p(p >= bounds(1) & p <= bounds(2));

    nexttile; hold on; box on; grid on

    histogram(p_trim,50,...
        'Normalization','probability',...
        'FaceColor',[0 0.447 0.741],...
        'EdgeColor','k');

    % True parameter value (red vertical line)
    xline(true_vals(i), 'r', 'LineWidth', 2);

    xlabel(['$',param_names{i},'$'],...
        'Interpreter','latex','FontSize',12);
    ylabel('Probability');
end