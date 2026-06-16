clc; 
clear;
close all;


data = readtable('Parameters_pop.xlsx', ...
    'ReadVariableNames', true, 'PreserveVariableNames', true);
data(1,:) = [];                          
num_individuals = height(data);


T0 = 167;  T_i0 = 85;  V0 = 69643;
y0 = [T0; T_i0; V0];
tspan = [0 500];


patientColor    = [0.00 0.35 0.95];      % blue for individuals
populationColor = [0.85 0.20 0.10];      % red for population
indLW = 0.5; popLW = 3.5;

%  Population parameters
lambda_pop  = 33.84;
k_pop       = 0.0087;
delta_1_pop = 0.00036;
delta_2_pop = 0.085;
p_pop       = 0.45;
c_pop       = 0.96;

% Axis styling helper 
set_axis_properties = @(ax) set(ax, ...
    'XTick',0:50:500, 'Box','on', 'LineWidth',1.0, 'FontSize',12, ...
    'XGrid','off', 'YGrid','on', 'GridLineStyle',':');

% Figure 1: Uninfected CD4+  T(t)
figure; hold on;
title('Uninfected CD4^{+} Cell Count Per Week');
xlabel('Time (weeks)', 'FontSize',12);
ylabel('$T(t)$', 'Interpreter','latex', 'FontSize',12);


for i = 1:num_individuals
    lambda  = data.lambda(i);  k = data.k(i);
    delta_1 = data.delta_1(i); delta_2 = data.delta_2(i);
    p       = data.p(i);       c = data.c(i);
    [t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda,k,delta_1,delta_2,p,c), tspan, y0);
    plot(t, y(:,1), 'Color', patientColor, 'LineWidth', indLW);
end
[t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda_pop,k_pop,delta_1_pop,delta_2_pop,p_pop,c_pop), tspan, y0);
plot(t, y(:,1), 'Color', populationColor, 'LineWidth', popLW);
set_axis_properties(gca);

%Figure 2: Infected CD4+  T_i(t)
figure; hold on;
title('Infected CD4^{+} Cell Counts Per Week');
xlabel('Time (weeks)', 'FontSize',12);
ylabel('$T_i(t)$', 'Interpreter','latex', 'FontSize',12);
ylim([0, 1800]);

for i = 1:num_individuals
    lambda  = data.lambda(i);  k = data.k(i);
    delta_1 = data.delta_1(i); delta_2 = data.delta_2(i);
    p       = data.p(i);       c = data.c(i);
    [t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda,k,delta_1,delta_2,p,c), tspan, y0);
    plot(t, y(:,2), 'Color', patientColor, 'LineWidth', indLW);
end
[t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda_pop,k_pop,delta_1_pop,delta_2_pop,p_pop,c_pop), tspan, y0);
plot(t, y(:,2), 'Color', populationColor, 'LineWidth', popLW);
set_axis_properties(gca);

% Figure 3: Total CD4+  T(t)+T_i(t)
figure; hold on;
title('Total CD4^{+} Cell Count Per Week');
xlabel('Time (weeks)', 'FontSize',12);
ylabel('$T(t)+T_i(t)$', 'Interpreter','latex', 'FontSize', 12);

for i = 1:num_individuals
    lambda  = data.lambda(i);  k = data.k(i);
    delta_1 = data.delta_1(i); delta_2 = data.delta_2(i);
    p       = data.p(i);       c = data.c(i);
    [t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda,k,delta_1,delta_2,p,c), tspan, y0);
    plot(t, y(:,1)+y(:,2), 'Color', patientColor, 'LineWidth', indLW);
end
[t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda_pop,k_pop,delta_1_pop,delta_2_pop,p_pop,c_pop), tspan, y0);
plot(t, y(:,1)+y(:,2), 'Color', populationColor, 'LineWidth', popLW);
set_axis_properties(gca);

% Figure 4: Viral load  V(t) 
figure; hold on;
title('Plasma HIV RNA Copies Per mL of Blood');
xlabel('Time (weeks)', 'FontSize',12);
ylabel('$V(t)$', 'Interpreter','latex', 'FontSize',12);

for i = 1:num_individuals
    lambda  = data.lambda(i);  k = data.k(i);
    delta_1 = data.delta_1(i); delta_2 = data.delta_2(i);
    p       = data.p(i);       c = data.c(i);
    [t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda,k,delta_1,delta_2,p,c), tspan, y0);
    plot(t, y(:,3), 'Color', patientColor, 'LineWidth', indLW);
end
[t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda_pop,k_pop,delta_1_pop,delta_2_pop,p_pop,c_pop), tspan, y0);
plot(t, y(:,3), 'Color', populationColor, 'LineWidth', popLW);

% Set log scale, turn OFF built-in grid
set(gca, 'YScale','log', ...
         'YTick',[10 100 1000 10000], ...  % ticks only at your chosen values
         'YGrid','off', 'XGrid','off', ...
         'YMinorGrid','off', 'YMinorTick','off'); % no minor lines at all

% Draw only the desired horizontal reference lines
arrayfun(@(v) yline(v,'--k','LineWidth',0.8), [10 100 1000 10000]);




[t, y] = ode23s(@(t,y) HIV_Within_Host(t,y,lambda_pop,k_pop,delta_1_pop,delta_2_pop,p_pop,c_pop), tspan, y0);
plot(t, y(:,3), 'Color', populationColor, 'LineWidth', popLW);
set(gca, 'YScale','log', 'YTick',[1e2 1e3 1e4], 'YMinorTick','off');
set_axis_properties(gca);

function dy = HIV_Within_Host(~, y, lambda, k, delta_1, delta_2, p, c)
    % y = [T; T_i; V]
    T = y(1); T_i = y(2); V = y(3);
    dy = zeros(3,1);
    dy(1) = lambda - k*T*V - delta_1*T;
    dy(2) = k*T*V - delta_2*T_i;
    dy(3) = p*T_i - c*V;
end
