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

stoch_simul(periods=200000,drop=2500, order=1,irf=100) y pi r epsilon_y epsilon_pi;

save illus_learning_sim y pi r epsilon_y epsilon_pi
