cd '/Users/julianashwin/Documents/GitHub/EcoNNet.jl/dynare_estimation/illustrative';

addpath /Applications/Dynare/6.0-x86_64/matlab
clear;
close all;

%% Set parameters 
% Determinate model 
beta = 0.95;
kappa = 0.05;
eta = 0.95;
sigma = 0.25;
phi_pi = 1.5;
sigma_pi = 0.2;
sigma_y = 0.2;
sigma_r = 0.01;
rho_pi = 0.5;
rho_y = 0.5;
save param_illus_det beta kappa eta sigma phi_pi sigma_pi sigma_y rho_pi rho_y
% Indeterminate model
beta = 0.95;
kappa = 0.05;
eta = 0.95;
sigma = 0.25;
phi_pi = 0.5;
sigma_pi = 0.2;
sigma_y = 0.2;
sigma_z = 0.2;
rho_pi = 0.5;
rho_y = 0.5;
omega12 = 0; % covariance of nu_y and nu_pi
omega13 = 0.3; % covariance of nu_y and sunspot
omega23 = 0.1; % covariance of nu_pi and sunspot

save param_illus_indet beta kappa eta sigma phi_pi  rho_pi rho_y ... 
    sigma_pi sigma_y sigma_z omega12 omega13 omega23


%% Save learning data

learning_data = readtable('illus_sim.csv');
y = learning_data.y;
pi = learning_data.pi;
p = learning_data.pi;
r = learning_data.r;
epsilon_y = learning_data.epsilon_y;
epsilon_pi = learning_data.epsilon_pi;

save illus_learning_sim y pi p r epsilon_y epsilon_pi



%% Determinate model

% Solve and simulate the linear determinate model
dynare illus_det noclearall nolog

% Estimate determinate model on data generated by dynare
dynare illus_det_est noclearall nolog

% Estimate on learning data
dynare illus_det_est_learn_martin noclearall nolog



%%  Set parameters and solve the indeterminate model with a sunspot shock

close all



% Solve and simulate the linear indeterminate model with sunspots
dynare illus_indet noclearall nolog

% Estimate the linear sunspot model on data generated by dynare
dynare illus_indet_est noclearall nolog

% Estimate the linear sunspot model on learning data
dynare illus_indet_est_learn noclearall nolog
