// Indeterminate NK model

var y pi ps epsilon_y epsilon_pi;
varexo nu_y nu_pi sunspot;

parameters beta kappa eta sigma phi_pi rho_pi rho_y
sigma_pi sigma_y sigma_z
omega12 omega13 omega23;


load param_illus_indet;
set_param_value('beta',beta); // discount
set_param_value('kappa',kappa); // output in PC
set_param_value('eta',eta); // persistence in y
set_param_value('sigma',sigma); // interest rate in y
set_param_value('phi_pi',phi_pi); // Taylor Rule
set_param_value('sigma_pi',sigma_pi); // inflation shock variance
set_param_value('sigma_y',sigma_y); // output shock variance
set_param_value('sigma_z',sigma_z); // ``sunspot" shock variance
set_param_value('omega12',omega12); // covariance of epsilon_y and epsilon_pi
set_param_value('omega13',omega13); // covariance of epsilon_y and sunspot
set_param_value('omega23',omega23); // covariance of epsilon_pi and sunspot
set_param_value('rho_pi',rho_pi); // inflation shock persistence
set_param_value('rho_y',rho_y); // output shock persistence




model(linear);
pi = beta*ps + kappa*y + epsilon_pi;
pi - ps(-1) = sunspot;
y = eta*y(-1) - sigma*(phi_pi*pi - ps) + epsilon_y;
epsilon_pi = rho_pi*epsilon_pi(-1) + nu_pi;
epsilon_y = rho_y*epsilon_y(-1) + nu_y;
end;


// Steady state
steady;

// Blanchard-Kahn conditions
check;

// Perturbation analysis
shocks;
var nu_pi; stderr sigma_pi;
var nu_y; stderr sigma_y;
var sunspot; stderr sigma_z;
corr nu_y, nu_pi = omega12;
corr nu_y, sunspot = omega13;
corr nu_pi, sunspot = omega23;
end;

varobs y pi;

estimated_params;
beta, 0.95,-3,3, normal_pdf, 1, 1,,;
kappa, 0.05,-3,3, normal_pdf, 0, 1,,;
eta, 0.95,-3,3, normal_pdf, 1, 1,,;
sigma, 0.25,-3,3, normal_pdf, 0, 1,,;
phi_pi, 0.5,0,5, normal_pdf, 2, 1,,;
rho_pi, 0,,, normal_pdf, 0, 0.5,,;
rho_y, 0,,, normal_pdf, 0, 0.5,,;

stderr nu_pi, 0.2, 0, 3, inv_gamma_pdf, 0.5, 2;
stderr nu_y, 0.2, 0, 3, inv_gamma_pdf, 0.5, 2;
stderr sunspot, 0.5, 0, 3, inv_gamma_pdf, 0.5, 2;
//corr nu_pi, nu_y,0,-1,1, normal_pdf, 0, 1;
corr nu_pi, sunspot, 0.5, -1, 1, beta_pdf, 0, 0.3, -1, 1;
corr nu_y, sunspot, 0.5, -1, 1, beta_pdf, 0, 0.3, -1, 1;

end;


estimation(datafile=illus_learning_temp,
first_obs = 1, nobs = 1000, mode_compute=6, mh_replic = 5000,
mh_jscale = 2, order = 1, nograph, nodisplay) y pi ;
