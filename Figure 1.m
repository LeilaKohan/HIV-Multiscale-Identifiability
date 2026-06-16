clc; clear; close all;

% Read the Monolix estimated individual parameters file
filename = 'simulatedIndividualParameters.txt';
opts = detectImportOptions(filename, 'Delimiter', ',');
opts.VariableNamesLine = 1;
opts.DataLine = 2;
data = readtable(filename, opts);

%data = readtable('Parameters_pop.xlsx', 'ReadVariableNames', true, 'PreserveVariableNames', true);
data(1, :) = [];  % Remove header




lambda_vals = data.lambda;
k_vals = data.beta;
delta_1_vals = data.d;
delta_2_vals = data.delta;
p_vals = data.pi;
c_vals = data.c;




param_names = {'\lambda', '\beta_v', '\delta_1', '\delta_2', 'p', 'c'};
param_vals = {lambda_vals, k_vals, delta_1_vals, delta_2_vals, p_vals, c_vals};


% Create figure
figure;
tiledlayout(2,3, 'TileSpacing', 'compact', 'Padding', 'compact');




for i = 1:length(param_vals)
    nexttile;
    vals = param_vals{i};
     
    histogram(vals, 'Normalization', 'probability', 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'k'); 
    hold on;
 
    xlabel(param_names{i}, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Probability', 'FontSize', 10);
    grid on;
end





