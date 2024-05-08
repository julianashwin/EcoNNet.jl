// Determinate two equation asset pricing model

var y p; 
varexo epsilon_d epsilon_s;

parameters phi_yy phi_yp phi_pp phi_py sigma_s sigma_d alpha_py;

load param_general_det;
set_param_value('phi_yy',phi_yy); // persistence of y_t
set_param_value('phi_yp',phi_yp); // feedback from p to y
set_param_value('phi_pp',phi_pp); // discount factor on future dividends
set_param_value('phi_py',phi_py); // feedback from output to dividends
set_param_value('sigma_s',sigma_s); // ``supply" shock variance
set_param_value('sigma_d',sigma_d); // ``demand" shock variance
set_param_value('alpha_py',alpha_py); // nonlinearity in p-y relationship
//set_param_value('rho_d',rho_d); // ``demand" shock persistence
//set_param_value('rho_s',rho_s); // ``supplu" shock persistence

model(linear); 
p = phi_pp*p(+1) + phi_py*y + epsilon_d;
y = phi_yy*y(-1) + phi_yp*p + epsilon_s;
//epsilon_d = rho_d*epsilon_d(-1) + nu_d;
//epsilon_s = rho_s*epsilon_s(-1) + nu_s;
end;

steady;
check;

shocks;
var epsilon_d; stderr sigma_d;
var epsilon_s; stderr sigma_s;
end;

varobs y p;

estimated_params;
phi_yy, 0.9,-2,2, normal_pdf, 1, 1,,;
phi_yp, 0.1,-2,2, normal_pdf, 0, 1,,;
phi_py, -0.85,-2,2, normal_pdf, 0, 1,,;
phi_pp, 0.95,-2,2, normal_pdf, 1, 1,,;
//rho_d, 0,-2,2; // normal_pdf, 0, 0.5,,;
//rho_s, 0,-2,2; // normal_pdf, 0, 0.5,,;

//stderr epsilon_d, 0.5, 0, 3, inv_gamma_pdf, 0.5, 2;
stderr epsilon_s, 0.5, 0, 3, inv_gamma_pdf, 0.5, 2;
//corr epsilon_s, epsilon_d,0,-1,1, normal_pdf, 0, 1;
end;

	
//set_dynare_seed('clock');

estimation(datafile=general_sim_det,
first_obs = 45000, mode_compute=2, mh_replic = 10000,
mh_jscale = 2, order = 1) y p ;


