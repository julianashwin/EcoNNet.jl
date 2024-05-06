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
stoch_simul(periods=50000,drop=2500, order=1,irf=20);


save 'general_sim_det' y p;