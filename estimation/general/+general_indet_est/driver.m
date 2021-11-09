%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'general_indet_est';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'epsilon_d'};
M_.exo_names_tex(1) = {'epsilon\_d'};
M_.exo_names_long(1) = {'epsilon_d'};
M_.exo_names(2) = {'epsilon_s'};
M_.exo_names_tex(2) = {'epsilon\_s'};
M_.exo_names_long(2) = {'epsilon_s'};
M_.exo_names(3) = {'sunspot'};
M_.exo_names_tex(3) = {'sunspot'};
M_.exo_names_long(3) = {'sunspot'};
M_.endo_names = cell(3,1);
M_.endo_names_tex = cell(3,1);
M_.endo_names_long = cell(3,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y'};
M_.endo_names_long(1) = {'y'};
M_.endo_names(2) = {'p'};
M_.endo_names_tex(2) = {'p'};
M_.endo_names_long(2) = {'p'};
M_.endo_names(3) = {'ps'};
M_.endo_names_tex(3) = {'ps'};
M_.endo_names_long(3) = {'ps'};
M_.endo_partitions = struct();
M_.param_names = cell(11,1);
M_.param_names_tex = cell(11,1);
M_.param_names_long = cell(11,1);
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
M_.param_names(5) = {'alpha_py'};
M_.param_names_tex(5) = {'alpha\_py'};
M_.param_names_long(5) = {'alpha_py'};
M_.param_names(6) = {'sigma_s'};
M_.param_names_tex(6) = {'sigma\_s'};
M_.param_names_long(6) = {'sigma_s'};
M_.param_names(7) = {'sigma_d'};
M_.param_names_tex(7) = {'sigma\_d'};
M_.param_names_long(7) = {'sigma_d'};
M_.param_names(8) = {'sigma_z'};
M_.param_names_tex(8) = {'sigma\_z'};
M_.param_names_long(8) = {'sigma_z'};
M_.param_names(9) = {'omega12'};
M_.param_names_tex(9) = {'omega12'};
M_.param_names_long(9) = {'omega12'};
M_.param_names(10) = {'omega13'};
M_.param_names_tex(10) = {'omega13'};
M_.param_names_long(10) = {'omega13'};
M_.param_names(11) = {'omega23'};
M_.param_names_tex(11) = {'omega23'};
M_.param_names_long(11) = {'omega23'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 3;
M_.param_nbr = 11;
M_.orig_endo_nbr = 3;
M_.aux_vars = [];
options_.varobs = cell(2, 1);
options_.varobs(1)  = {'y'};
options_.varobs(2)  = {'p'};
options_.varobs_id = [ 1 2  ];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
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
M_.orig_eq_nbr = 3;
M_.eq_nbr = 3;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 0;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 0;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 3;
 0 4;
 2 5;]';
M_.nstatic = 1;
M_.nfwrd   = 0;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 0;
M_.nspred   = 2;
M_.ndynamic   = 2;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'p' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'y' ;
};
M_.mapping.y.eqidx = [1 3 ];
M_.mapping.p.eqidx = [1 2 3 ];
M_.mapping.ps.eqidx = [1 2 ];
M_.mapping.epsilon_d.eqidx = [1 ];
M_.mapping.epsilon_s.eqidx = [3 ];
M_.mapping.sunspot.eqidx = [2 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 3 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 0;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 0;
oo_.steady_state = zeros(3, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(11, 1);
M_.endo_trends = struct('deflator', cell(3, 1), 'log_deflator', cell(3, 1), 'growth_factor', cell(3, 1), 'log_growth_factor', cell(3, 1));
M_.NNZDerivatives = [11; 0; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load param_general_indet;
set_param_value('phi_yy',phi_yy); 
set_param_value('phi_yp',phi_yp); 
set_param_value('phi_pp',phi_pp); 
set_param_value('phi_py',phi_py); 
set_param_value('sigma_s',sigma_s); 
set_param_value('sigma_d',sigma_d); 
set_param_value('sigma_z',sigma_z); 
set_param_value('alpha_py',alpha_py); 
set_param_value('omega12',omega12); 
set_param_value('omega13',omega13); 
set_param_value('omega23',omega23); 
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(7))^2;
M_.Sigma_e(2, 2) = (M_.params(6))^2;
M_.Sigma_e(3, 3) = (M_.params(8))^2;
M_.Sigma_e(1, 3) = M_.params(11)*sqrt(M_.Sigma_e(1, 1)*M_.Sigma_e(3, 3));
M_.Sigma_e(3, 1) = M_.Sigma_e(1, 3);
M_.Correlation_matrix(1, 3) = M_.params(11);
M_.Correlation_matrix(3, 1) = M_.Correlation_matrix(1, 3);
M_.Sigma_e(2, 1) = M_.params(9)*sqrt(M_.Sigma_e(2, 2)*M_.Sigma_e(1, 1));
M_.Sigma_e(1, 2) = M_.Sigma_e(2, 1);
M_.Correlation_matrix(2, 1) = M_.params(9);
M_.Correlation_matrix(1, 2) = M_.Correlation_matrix(2, 1);
M_.Sigma_e(2, 3) = M_.params(10)*sqrt(M_.Sigma_e(2, 2)*M_.Sigma_e(3, 3));
M_.Sigma_e(3, 2) = M_.Sigma_e(2, 3);
M_.Correlation_matrix(2, 3) = M_.params(10);
M_.Correlation_matrix(3, 2) = M_.Correlation_matrix(2, 3);
M_.sigma_e_is_diagonal = 0;
estim_params_.var_exo = zeros(0, 10);
estim_params_.var_endo = zeros(0, 10);
estim_params_.corrx = zeros(0, 11);
estim_params_.corrn = zeros(0, 11);
estim_params_.param_vals = zeros(0, 10);
estim_params_.param_vals = [estim_params_.param_vals; 1, 0.5, (-2), 2, 3, 1, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 2, 0.4, (-2), 2, 3, 0, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, 0.5, (-2), 2, 3, 0, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 3, 0.95, (-2), 2, 3, 1, 1, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 2, 0.5, 0, 3, 4, 0.5, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 3, 0.5, 0, 3, 4, 0.5, 2, NaN, NaN, NaN ];
estim_params_.corrx = [estim_params_.corrx; 1, 3, 0.5, (-1), 1, 1, 0, 0.3, (-1), 1, NaN ];
estim_params_.corrx = [estim_params_.corrx; 2, 3, 0.5, (-1), 1, 1, 0, 0.3, (-1), 1, NaN ];
options_.mh_jscale = 2;
options_.mh_replic = 10000;
options_.mode_compute = 2;
options_.order = 1;
options_.datafile = 'general_sim_indet';
options_.first_obs = 45000;
var_list_ = {'y';'p'};
oo_recursive_=dynare_estimation(var_list_);
save('general_indet_est_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('general_indet_est_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('general_indet_est_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('general_indet_est_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('general_indet_est_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('general_indet_est_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('general_indet_est_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
