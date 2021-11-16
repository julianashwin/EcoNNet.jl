// Indeterminate two equation asset pricing model

var y p ps epsilon_d epsilon_s ; 
varexo nu_d nu_s sunspot;

parameters phi_yy phi_yp phi_pp phi_py alpha_py rho_d rho_s
sigma_s sigma_d sigma_z
omega12 omega13 omega23;

load param_general_indet;
set_param_value('phi_yy',phi_yy); // persistence of y_t
set_param_value('phi_yp',phi_yp); // feedback from p to y
set_param_value('phi_pp',phi_pp); // discount factor on future dividends
set_param_value('phi_py',phi_py); // feedback from output to dividends
set_param_value('rho_s',rho_s); // ``supply" shock persistence
set_param_value('rho_d',rho_d); // ``demand" shock persistence
set_param_value('sigma_s',sigma_s); // ``supply" shock variance
set_param_value('sigma_d',sigma_d); // ``demand" shock variance
set_param_value('sigma_z',sigma_z); // ``sunspot" shock variance
set_param_value('alpha_py',alpha_py); // nonlinearity in p-y relationship
set_param_value('omega12',omega12); // covariance of epsilon_s and epsilon_d
set_param_value('omega13',omega13); // covariance of epsilon_s and sunspot
set_param_value('omega23',omega23); // covariance of epsilon_d and sunspot


model(linear); 
p = phi_pp*ps + phi_py*y + epsilon_d;
p - ps(-1) = sunspot;
y = phi_yy*y(-1) + phi_yp*p + epsilon_s;
epsilon_s = rho_s*epsilon_s(-1) + nu_s;
epsilon_d = rho_d*epsilon_d(-1) + nu_d;
end;

// Steady state
steady;

// Blanchard-Kahn conditions
check;

// Perturbation analysis
shocks;

var nu_s; stderr sigma_s;
var nu_d; stderr sigma_d;
var sunspot; stderr sigma_z;
corr nu_s, nu_d = omega12;
corr nu_s, sunspot = omega13;
corr nu_d, sunspot = omega23;
end;

varobs y p;

estimated_params;
phi_yy, 0.5,-2,2, normal_pdf, 1, 1,,;
phi_yp, 0.4,-2,2, normal_pdf, 0, 1,,;
phi_py, 0.5,-2,2, normal_pdf, 0, 1,,;
phi_pp, 0.95,-2,2, normal_pdf, 1, 1,,;
rho_d, 0,-2,2, normal_pdf, 0, 0.5,,;
rho_s, 0,-2,2, normal_pdf, 0, 0.5,,;

stderr nu_d, 0.5, 0, 3, inv_gamma_pdf, 0.5, 2;
stderr nu_s, 0.5, 0, 3, inv_gamma_pdf, 0.5, 2;
stderr sunspot, 0.5, 0, 3, inv_gamma_pdf, 0.5, 2;
//corr nu_s, nu_d,0,-1,1, normal_pdf, 0, 1;
corr nu_d, sunspot, 0.5, -1, 1, beta_pdf, 0, 0.3, -1, 1;
corr nu_s, sunspot, 0.5, -1, 1, beta_pdf, 0, 0.3, -1, 1;

end;


estimation(datafile=illus_learning_sim,
first_obs = 90000, mode_compute=6, mh_replic = 5000,
mh_jscale = 2, order = 1) y p;
