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
stoch_simul(periods=100000,drop=2500, order=1,irf=100);


save 'illus_sim_indet' y pi;

