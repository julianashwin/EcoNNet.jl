%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'illus_indet';
M_.dynare_version = '6.0';
oo_.dynare_version = '6.0';
options_.dynare_version = '6.0';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'nu_y'};
M_.exo_names_tex(1) = {'nu\_y'};
M_.exo_names_long(1) = {'nu_y'};
M_.exo_names(2) = {'nu_pi'};
M_.exo_names_tex(2) = {'nu\_pi'};
M_.exo_names_long(2) = {'nu_pi'};
M_.exo_names(3) = {'sunspot'};
M_.exo_names_tex(3) = {'sunspot'};
M_.exo_names_long(3) = {'sunspot'};
M_.endo_names = cell(5,1);
M_.endo_names_tex = cell(5,1);
M_.endo_names_long = cell(5,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y'};
M_.endo_names_long(1) = {'y'};
M_.endo_names(2) = {'pi'};
M_.endo_names_tex(2) = {'pi'};
M_.endo_names_long(2) = {'pi'};
M_.endo_names(3) = {'ps'};
M_.endo_names_tex(3) = {'ps'};
M_.endo_names_long(3) = {'ps'};
M_.endo_names(4) = {'epsilon_y'};
M_.endo_names_tex(4) = {'epsilon\_y'};
M_.endo_names_long(4) = {'epsilon_y'};
M_.endo_names(5) = {'epsilon_pi'};
M_.endo_names_tex(5) = {'epsilon\_pi'};
M_.endo_names_long(5) = {'epsilon_pi'};
M_.endo_partitions = struct();
M_.param_names = cell(13,1);
M_.param_names_tex = cell(13,1);
M_.param_names_long = cell(13,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'kappa'};
M_.param_names_tex(2) = {'kappa'};
M_.param_names_long(2) = {'kappa'};
M_.param_names(3) = {'eta'};
M_.param_names_tex(3) = {'eta'};
M_.param_names_long(3) = {'eta'};
M_.param_names(4) = {'sigma'};
M_.param_names_tex(4) = {'sigma'};
M_.param_names_long(4) = {'sigma'};
M_.param_names(5) = {'phi_pi'};
M_.param_names_tex(5) = {'phi\_pi'};
M_.param_names_long(5) = {'phi_pi'};
M_.param_names(6) = {'rho_pi'};
M_.param_names_tex(6) = {'rho\_pi'};
M_.param_names_long(6) = {'rho_pi'};
M_.param_names(7) = {'rho_y'};
M_.param_names_tex(7) = {'rho\_y'};
M_.param_names_long(7) = {'rho_y'};
M_.param_names(8) = {'sigma_pi'};
M_.param_names_tex(8) = {'sigma\_pi'};
M_.param_names_long(8) = {'sigma_pi'};
M_.param_names(9) = {'sigma_y'};
M_.param_names_tex(9) = {'sigma\_y'};
M_.param_names_long(9) = {'sigma_y'};
M_.param_names(10) = {'sigma_z'};
M_.param_names_tex(10) = {'sigma\_z'};
M_.param_names_long(10) = {'sigma_z'};
M_.param_names(11) = {'omega12'};
M_.param_names_tex(11) = {'omega12'};
M_.param_names_long(11) = {'omega12'};
M_.param_names(12) = {'omega13'};
M_.param_names_tex(12) = {'omega13'};
M_.param_names_long(12) = {'omega13'};
M_.param_names(13) = {'omega23'};
M_.param_names_tex(13) = {'omega23'};
M_.param_names_long(13) = {'omega23'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 5;
M_.param_nbr = 13;
M_.orig_endo_nbr = 5;
M_.aux_vars = [];
options_.varobs = cell(2, 1);
options_.varobs(1)  = {'y'};
options_.varobs(2)  = {'pi'};
options_.varobs_id = [ 1 2  ];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.learnt_shocks = [];
M_.learnt_endval = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
M_.matched_irfs = {};
M_.matched_irfs_weights = {};
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.eq_nbr = 5;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
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
 1 5;
 0 6;
 2 7;
 3 8;
 4 9;]';
M_.nstatic = 1;
M_.nfwrd   = 0;
M_.npred   = 4;
M_.nboth   = 0;
M_.nsfwrd   = 0;
M_.nspred   = 4;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , 'pi' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'y' ;
  4 , 'name' , 'epsilon_pi' ;
  5 , 'name' , 'epsilon_y' ;
};
M_.mapping.y.eqidx = [1 3 ];
M_.mapping.pi.eqidx = [1 2 3 ];
M_.mapping.ps.eqidx = [1 2 3 ];
M_.mapping.epsilon_y.eqidx = [3 5 ];
M_.mapping.epsilon_pi.eqidx = [1 4 ];
M_.mapping.nu_y.eqidx = [5 ];
M_.mapping.nu_pi.eqidx = [4 ];
M_.mapping.sunspot.eqidx = [2 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = true;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 3;
M_.block_structure.block(1).mfs = 3;
M_.block_structure.block(1).equation = [ 2 4 5];
M_.block_structure.block(1).variable = [ 2 5 4];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 3;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [2 3 4 5 6 ];
M_.block_structure.block(2).Simulation_Type = 6;
M_.block_structure.block(2).endo_nbr = 2;
M_.block_structure.block(2).mfs = 1;
M_.block_structure.block(2).equation = [ 1 3];
M_.block_structure.block(2).variable = [ 3 1];
M_.block_structure.block(2).is_linear = true;
M_.block_structure.block(2).NNZDerivatives = 2;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [0 0 1 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([1 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([1 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure.variable_reordered = [ 2 5 4 3 1];
M_.block_structure.equation_reordered = [ 2 4 5 1 3];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 2 3;
 3 1;
 4 5;
 5 4;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 1;
 1 2;
 1 3;
 1 5;
 2 2;
 3 1;
 3 2;
 3 3;
 3 4;
 4 5;
 5 4;
];
M_.block_structure.dyn_tmp_nbr = 0;
M_.state_var = [5 4 3 1 ];
M_.maximum_lag = 1;
M_.maximum_lead = 0;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 0;
oo_.steady_state = zeros(5, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(13, 1);
M_.endo_trends = struct('deflator', cell(5, 1), 'log_deflator', cell(5, 1), 'growth_factor', cell(5, 1), 'log_growth_factor', cell(5, 1));
M_.NNZDerivatives = [18; 0; -1; ];
M_.dynamic_g1_sparse_rowval = int32([3 2 5 4 1 3 1 2 3 1 3 3 5 1 4 5 4 2 ]);
M_.dynamic_g1_sparse_colval = int32([1 3 4 5 6 6 7 7 7 8 8 9 9 10 10 16 17 18 ]);
M_.dynamic_g1_sparse_colptr = int32([1 2 2 3 4 5 7 10 12 14 16 16 16 16 16 16 17 18 19 ]);
M_.dynamic_g2_sparse_indices = int32([]);
M_.lhs = {
'pi'; 
'pi-ps(-1)'; 
'y'; 
'epsilon_pi'; 
'epsilon_y'; 
};
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 4];
M_.block_structure_stat.block(1).variable = [ 5];
M_.block_structure_stat.block(2).Simulation_Type = 3;
M_.block_structure_stat.block(2).endo_nbr = 1;
M_.block_structure_stat.block(2).mfs = 1;
M_.block_structure_stat.block(2).equation = [ 5];
M_.block_structure_stat.block(2).variable = [ 4];
M_.block_structure_stat.block(3).Simulation_Type = 6;
M_.block_structure_stat.block(3).endo_nbr = 3;
M_.block_structure_stat.block(3).mfs = 3;
M_.block_structure_stat.block(3).equation = [ 3 1 2];
M_.block_structure_stat.block(3).variable = [ 1 2 3];
M_.block_structure_stat.variable_reordered = [ 5 4 1 2 3];
M_.block_structure_stat.equation_reordered = [ 4 5 3 1 2];
M_.block_structure_stat.incidence.sparse_IM = [
 1 1;
 1 2;
 1 3;
 1 5;
 2 2;
 2 3;
 3 1;
 3 2;
 3 3;
 3 4;
 4 5;
 5 4;
];
M_.block_structure_stat.tmp_nbr = 0;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([1 2 1 2 3 1 2 3 ]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([1 1 2 2 2 3 3 3 ]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([1 3 6 9 ]);
M_.static_g1_sparse_rowval = int32([1 3 1 2 3 1 2 3 3 5 1 4 ]);
M_.static_g1_sparse_colval = int32([1 1 2 2 2 3 3 3 4 4 5 5 ]);
M_.static_g1_sparse_colptr = int32([1 3 6 9 11 13 ]);
load param_illus_indet;
set_param_value('beta',beta); 
set_param_value('kappa',kappa); 
set_param_value('eta',eta); 
set_param_value('sigma',sigma); 
set_param_value('phi_pi',phi_pi); 
set_param_value('sigma_pi',sigma_pi); 
set_param_value('sigma_y',sigma_y); 
set_param_value('sigma_z',sigma_z); 
set_param_value('omega12',omega12); 
set_param_value('omega13',omega13); 
set_param_value('omega23',omega23); 
set_param_value('rho_pi',rho_pi); 
set_param_value('rho_y',rho_y); 
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(9))^2;
M_.Sigma_e(2, 2) = (M_.params(8))^2;
M_.Sigma_e(3, 3) = (M_.params(10))^2;
M_.Sigma_e(1, 2) = M_.params(11)*sqrt(M_.Sigma_e(1, 1)*M_.Sigma_e(2, 2));
M_.Sigma_e(2, 1) = M_.Sigma_e(1, 2);
M_.Correlation_matrix(1, 2) = M_.params(11);
M_.Correlation_matrix(2, 1) = M_.Correlation_matrix(1, 2);
M_.Sigma_e(1, 3) = M_.params(12)*sqrt(M_.Sigma_e(1, 1)*M_.Sigma_e(3, 3));
M_.Sigma_e(3, 1) = M_.Sigma_e(1, 3);
M_.Correlation_matrix(1, 3) = M_.params(12);
M_.Correlation_matrix(3, 1) = M_.Correlation_matrix(1, 3);
M_.Sigma_e(2, 3) = M_.params(13)*sqrt(M_.Sigma_e(2, 2)*M_.Sigma_e(3, 3));
M_.Sigma_e(3, 2) = M_.Sigma_e(2, 3);
M_.Correlation_matrix(2, 3) = M_.params(13);
M_.Correlation_matrix(3, 2) = M_.Correlation_matrix(2, 3);
M_.sigma_e_is_diagonal = 0;
options_.drop = 2500;
options_.irf = 100;
options_.order = 1;
options_.periods = 100000;
options_.replic = 50;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save 'illus_sim_indet' y pi;


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_indet_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
