cd '/Users/julianashwin/Documents/GitHub/EcoNNet.jl/estimation/illustrative';

addpath /Applications/Dynare/4.6.4/matlab
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

beta = -0.0621;
kappa = 0.8184;
eta = 0.9673;
sigma = -0.0154;
phi_pi = 2.1320;
sigma_pi = 0.4707;
sigma_y = 0.2157;
rho_pi = 0.5880;
rho_y = 0.4707;

save param_illus_det beta kappa eta sigma phi_pi sigma_pi sigma_y rho_pi rho_y
% Indeterminate model
beta = 1.06;
kappa = -0.03;
eta = 0.85;
sigma = 0.41;
phi_pi = 0.55;
sigma_pi = 0.25;
sigma_y = 0.17;
sigma_z = 0.39;
rho_pi = 0.45;
rho_y = 0.50;
omega12 = 0; % covariance of nu_y and nu_pi
omega13 = 0.70; % covariance of nu_y and sunspot
omega23 = 0.51; % covariance of nu_pi and sunspot

save param_illus_indet beta kappa eta sigma phi_pi  rho_pi rho_y ... 
    sigma_pi sigma_y sigma_z omega12 omega13 omega23




%% Solve models
close all 
% Solve and simulate the linear determinate model
dynare illus_det noclearall nolog

% Solve and simulate the linear indeterminate model with sunspots
dynare illus_indet noclearall nolog

