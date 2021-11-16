%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'general_det_est_learn';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'nu_d'};
M_.exo_names_tex(1) = {'nu\_d'};
M_.exo_names_long(1) = {'nu_d'};
M_.exo_names(2) = {'nu_s'};
M_.exo_names_tex(2) = {'nu\_s'};
M_.exo_names_long(2) = {'nu_s'};
M_.endo_names = cell(4,1);
M_.endo_names_tex = cell(4,1);
M_.endo_names_long = cell(4,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y'};
M_.endo_names_long(1) = {'y'};
M_.endo_names(2) = {'p'};
M_.endo_names_tex(2) = {'p'};
M_.endo_names_long(2) = {'p'};
M_.endo_names(3) = {'epsilon_d'};
M_.endo_names_tex(3) = {'epsilon\_d'};
M_.endo_names_long(3) = {'epsilon_d'};
M_.endo_names(4) = {'epsilon_s'};
M_.endo_names_tex(4) = {'epsilon\_s'};
M_.endo_names_long(4) = {'epsilon_s'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
M_.param_names(1) = {'phi_yy'};
M_.param_names_tex(1) = {'phi\_yy'};
M_.param_names_long(1) = {'phi_yy'};
M_.param_names(2) = {'phi_yp'};
M_.param_names_tex(2) = {'phi\_yp'};
M_.param_names_long(2) = {'phi_yp'};
M_.param_names(3) = {'phi_pp'};
M_.param_names_tex(3) = {'phi\_pp'};
M_.param_names_long(3) = {'phi_pp'};
M_.param_names(4) = {'phi_py'};
M_.param_names_tex(4) = {'phi\_py'};
M_.param_names_long(4) = {'phi_py'};
M_.param_names(5) = {'sigma_s'};
M_.param_names_tex(5) = {'sigma\_s'};
M_.param_names_long(5) = {'sigma_s'};
M_.param_names(6) = {'sigma_d'};
M_.param_names_tex(6) = {'sigma\_d'};
M_.param_names_long(6) = {'sigma_d'};
M_.param_names(7) = {'alpha_py'};
M_.param_names_tex(7) = {'alpha\_py'};
M_.param_names_long(7) = {'alpha_py'};
M_.param_names(8) = {'rho_d'};
M_.param_names_tex(8) = {'rho\_d'};
M_.param_names_long(8) = {'rho_d'};
M_.param_names(9) = {'rho_s'};
M_.param_names_tex(9) = {'rho\_s'};
M_.param_names_long(9) = {'rho_s'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 4;
M_.param_nbr = 9;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
options_.varobs = cell(2, 1);
options_.varobs(1)  = {'y'};
options_.varobs(2)  = {'p'};
options_.varobs_id = [ 1 2  ];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 4;
M_.eq_nbr = 4;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 4 0;
 0 5 8;
 2 6 0;
 3 7 0;]';
M_.nstatic = 0;
M_.nfwrd   = 1;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 1;
M_.nspred   = 3;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'p' ;
  2 , 'name' , 'y' ;
  3 , 'name' , 'epsilon_d' ;
  4 , 'name' , 'epsilon_s' ;
};
M_.mapping.y.eqidx = [1 2 ];
M_.mapping.p.eqidx = [1 2 ];
M_.mapping.epsilon_d.eqidx = [1 3 ];
M_.mapping.epsilon_s.eqidx = [2 4 ];
M_.mapping.nu_d.eqidx = [3 ];
M_.mapping.nu_s.eqidx = [4 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 3 4 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(9, 1);
M_.endo_trends = struct('deflator', cell(4, 1), 'log_deflator', cell(4, 1), 'growth_factor', cell(4, 1), 'log_growth_factor', cell(4, 1));
M_.NNZDerivatives = [14; 0; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load param_general_det;
set_param_value('phi_yy',phi_yy); 
set_param_value('phi_yp',phi_yp); 
set_param_value('phi_pp',phi_pp); 
set_param_value('phi_py',phi_py); 
set_param_value('sigma_s',sigma_s); 
set_param_value('sigma_d',sigma_d); 
set_param_value('alpha_py',alpha_py); 
set_param_value('rho_d',rho_d); 
set_param_value('rho_s',rho_s); 
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(6))^2;
M_.Sigma_e(2, 2) = (M_.params(5))^2;
estim_params_.var_exo = zeros(0, 10);
estim_params_.var_endo = zeros(0, 10);
estim_params_.corrx = zeros(0, 11);
estim_params_.corrn = zeros(0, 11);
estim_params_.param_vals = zeros(0, 10);
estim_params_.param_vals = [estim_params_.param_vals; 1, 0.9, (-2), 2, 3, 1, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 2, 0.1, (-2), 2, 3, 0, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, (-0.85), (-2), 2, 3, 0, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 3, 0.95, (-2), 2, 3, 1, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 8, 0, (-2), 2, 3, 0, 0.5, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 9, 0, (-2), 2, 3, 0, 0.5, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 1, 0.5, 0, 3, 4, 0.5, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 2, 0.5, 0, 3, 4, 0.5, 2, NaN, NaN, NaN ];
options_.mh_jscale = 2;
options_.mh_replic = 5000;
options_.mode_compute = 6;
options_.order = 1;
options_.datafile = 'illus_learning_sim';
options_.first_obs = 90000;
var_list_ = {'y';'p'};
oo_recursive_=dynare_estimation(var_list_);
save('general_det_est_learn_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('general_det_est_learn_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('general_det_est_learn_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('general_det_est_learn_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('general_det_est_learn_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('general_det_est_learn_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('general_det_est_learn_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
