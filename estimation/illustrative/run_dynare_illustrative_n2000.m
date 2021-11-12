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


%% Loop through many samples

learning_data = readtable('illus_sim.csv');

nruns = 10;
samplesize = 2000;
Run = transpose(1:nruns);
StartObs = transpose(1:samplesize:(nruns*samplesize));
SampleSize = ones(nruns,1)*samplesize;
DetDensity = zeros(nruns,1);
IndetDensity = zeros(nruns,1);
results_tab = table(Run,StartObs,SampleSize,DetDensity,IndetDensity);

for ii = 1:nruns
    display("Running estimation "+string(ii))
    close all
    % Select relevant data
    startob = results_tab.StartObs(ii);
    ssize = results_tab.SampleSize(ii); 
    y = learning_data.y(startob:(startob+ssize));
    pi = learning_data.pi(startob:(startob+ssize));
    r = learning_data.r(startob:(startob+ssize));
    epsilon_y = learning_data.epsilon_y(startob:(startob+ssize));
    epsilon_pi = learning_data.epsilon_pi(startob:(startob+ssize));
    % Save this slice of the data
    save illus_learning_temp y pi r epsilon_y epsilon_pi
    try
        % Estimate determinate model
        dynare illus_det_est_learn_n2000 noclearall nolog
        results_tab.DetDensity(ii) = oo_.MarginalDensity.ModifiedHarmonicMean;
    
        % Estimate indeterminate model
        dynare illus_indet_est_learn_n2000 noclearall nolog
        results_tab.IndetDensity(ii) = oo_.MarginalDensity.ModifiedHarmonicMean;
    end
    % Save results to csv
    writetable(results_tab,'results_table_'+string(ssize)+'.csv','Delimiter',',','QuoteStrings',true)
    endrred
