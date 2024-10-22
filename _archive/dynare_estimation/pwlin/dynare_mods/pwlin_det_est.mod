// Determinate NK model

var y pi epsilon_y epsilon_pi;
varexo nu_y nu_pi;

parameters beta kappa eta sigma phi_pi sigma_pi sigma_y rho_pi rho_y;

load param_illus_det;
set_param_value('beta',beta); // discount
set_param_value('kappa',kappa); // output in PC
set_param_value('eta',eta); // persistence in y
set_param_value('sigma',sigma); // interest rate in y
set_param_value('phi_pi',phi_pi); // Taylor Rule
set_param_value('sigma_pi',sigma_pi); // inflation shock variance
set_param_value('sigma_y',sigma_y); // output shock variance
set_param_value('rho_pi',rho_pi); // inflation shock persistence
set_param_value('rho_y',rho_y); // output shock persistence

model(linear);
pi = beta*pi(+1) + kappa*y + epsilon_pi;
y = eta*y(-1) - sigma*(phi_pi*pi - pi(+1)) + epsilon_y;
epsilon_pi = rho_pi*epsilon_pi(-1) + nu_pi;
epsilon_y = rho_y*epsilon_y(-1) + nu_y;
end;

steady;
check;

shocks;
var nu_pi; stderr sigma_pi;
var nu_y; stderr sigma_y;
end;

varobs y pi;

estimated_params;
beta, 0.95,-2,2, normal_pdf, 1, 1,,;
kappa, 0.05,-2,2, normal_pdf, 0, 1,,;
eta, 0.95,-2,2, normal_pdf, 1, 1,,;
sigma, 0.25,-2,2, normal_pdf, 0, 1,,;
phi_pi, 1.5,-2,2, normal_pdf, 2, 1,,;
//rho_d, 0,-2,2; // normal_pdf, 0, 0.5,,;
//rho_s, 0,-2,2; // normal_pdf, 0, 0.5,,;

//stderr nu_pi, 0.2, 0, 3, inv_gamma_pdf, 0.5, 2;
stderr nu_y, 0.2, 0, 3, inv_gamma_pdf, 0.5, 2;
//corr nu_pi, nu_y,0,-1,1, normal_pdf, 0, 1;
end;

	
//set_dynare_seed('clock');

estimation(datafile=illus_sim_det,
first_obs = 95000, mode_compute=6, mh_replic = 10000,
mh_jscale = 2, order = 1) y pi ;


