clear all
close all
clc

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

  %params = [48.9021714511231, 0.0176, 0.145294069286879, 0.2, 0.00185805364116262];
 params = [8.38663746616526	0.306293687587148	0.120967014736454	0.200000005673198	0.00119503317441228];
 nonlcon = @constraints_betas;
 k = fmincon(@(k) err_in_data(k), params, [], [], [], [], lb, ub, nonlcon, optimset('Display','iter'));

% results 

fprintf('Lambda = %g\n',  k(1));   
fprintf('beta = %g\n', k(2));
fprintf('gamma = %g\n',  k(3));
fprintf('rho = %g\n',  k(4));
%fprintf('alpha = %g\n',  k(5));
fprintf('beta_d = %g\n', k(5));

[HIV_Cases, HIV_Diagnose, AIDS_Cases, AIDS_deaths, I, D] = FiniteDifference_HIV_Between_Host_Model(k);

base_year = 1981;
real_tforward = tforward + base_year;
real_tforward_trimmed = tforward(1:end-1) + base_year;
real_year_of_infection = year_of_infection + base_year;
real_years_diagnose = years_diagnose + base_year;
real_years_AIDS_classification = years_AIDS_classification + base_year;
real_years_AIDS_death = years_AIDS_death + base_year;
% (a) Estimated HIV Incidence Over Time
figure;
plot(real_tforward_trimmed, HIV_Cases, '-b', 'LineWidth', 2); hold on;
plot(real_year_of_infection, Estimated_HIV_incident_numbers, 'r.', 'MarkerSize', 25);
xlabel('Year');
ylabel('HIV Incidence (in thousands)');
title('(a) Estimated Annual HIV Incidence in the U.S.');
%legend('Model', 'CDC Data', 'Location', 'NorthEast');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);  

% (b) Diagnosed HIV Population Over Time
figure;
plot(real_tforward, HIV_Diagnose, '-b', 'LineWidth', 2); hold on;
plot(real_years_diagnose, HIV_diagnose_cases, 'r.', 'MarkerSize', 25);
xlabel('Year');
ylabel('HIV Diagnoses (in thousands)');
title('(b) Annual HIV Diagnoses in the U.S.');
%legend('Model', 'CDC Data', 'Location', 'NorthEast');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);  

% (c) AIDS Classifications Over Time
figure;
plot(real_tforward_trimmed, AIDS_Cases, '-b', 'LineWidth', 2); hold on;
plot(real_years_AIDS_classification, AIDS_classidication_cases, 'r.', 'MarkerSize', 25);
xlabel('Year');
ylabel('AIDS Classifications (in thousands)');
title('(c) Annual AIDS Classifications in the U.S.');
%legend('Model', 'CDC Data', 'Location', 'NorthEast');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);  



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

function error_in_data = err_in_data(k)



Estimated_HIV_incident_numbers = [20.77, 65.39, 65.13, 131.07, 131.37,...
    85.37, 85.15, 49.53, 49.89, 58.6, 59.19, 56.40, 56.49, 53.70,...
    42.21, 41.10, 40.28, 38.89, 38.07, 38.66, 38.95,...
    38.68, 37.85, 37.30, 36.19]'; % Number of HIV incidents in 1000s



HIV_diagnose_cases = [46703, 44324, 42775, 40945, 40183, 38931, 39653,...
    39543, 39163, 38048, 37132, 36350]'/1000; % CDC DATA




AIDS_classidication_cases = [11808, 18806, 27929, 34432, 41588, 47387, 57222,...
    69709, 73971, 67469, 64353, 56500, 45946, 39343, 37728, 38049, 36622,...
    36429, 36992, 35900, 33956, 32616, 31821, 31234, 30063, 27328, 25556,...
    24608, 23574, 19244, 18539, 18284, 17727, 17106, 16410]'/1000; % CDC 




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



error_in_data = 2*sum((Model_HIV_Inc -  Estimated_HIV_incident_numbers).^2)+...
                20*sum((Model_HIV_Diagnoses - HIV_diagnose_cases).^2) +...
                11*sum((Model_AIDS_Cases - AIDS_classidication_cases).^2);
              


              
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
