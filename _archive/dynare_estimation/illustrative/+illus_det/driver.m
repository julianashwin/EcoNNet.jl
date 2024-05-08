%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'illus_det';
M_.dynare_version = '6.0';
oo_.dynare_version = '6.0';
options_.dynare_version = '6.0';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'nu_y'};
M_.exo_names_tex(1) = {'nu\_y'};
M_.exo_names_long(1) = {'nu_y'};
M_.exo_names(2) = {'nu_pi'};
M_.exo_names_tex(2) = {'nu\_pi'};
M_.exo_names_long(2) = {'nu_pi'};
M_.endo_names = cell(4,1);
M_.endo_names_tex = cell(4,1);
M_.endo_names_long = cell(4,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y'};
M_.endo_names_long(1) = {'y'};
M_.endo_names(2) = {'pi'};
M_.endo_names_tex(2) = {'pi'};
M_.endo_names_long(2) = {'pi'};
M_.endo_names(3) = {'epsilon_y'};
M_.endo_names_tex(3) = {'epsilon\_y'};
M_.endo_names_long(3) = {'epsilon_y'};
M_.endo_names(4) = {'epsilon_pi'};
M_.endo_names_tex(4) = {'epsilon\_pi'};
M_.endo_names_long(4) = {'epsilon_pi'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
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
M_.param_names(6) = {'sigma_pi'};
M_.param_names_tex(6) = {'sigma\_pi'};
M_.param_names_long(6) = {'sigma_pi'};
M_.param_names(7) = {'sigma_y'};
M_.param_names_tex(7) = {'sigma\_y'};
M_.param_names_long(7) = {'sigma_y'};
M_.param_names(8) = {'rho_pi'};
M_.param_names_tex(8) = {'rho\_pi'};
M_.param_names_long(8) = {'rho_pi'};
M_.param_names(9) = {'rho_y'};
M_.param_names_tex(9) = {'rho\_y'};
M_.param_names_long(9) = {'rho_y'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 4;
M_.param_nbr = 9;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
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
M_.eq_nbr = 4;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
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
M_.equations_tags = {
  1 , 'name' , 'pi' ;
  2 , 'name' , 'y' ;
  3 , 'name' , 'epsilon_pi' ;
  4 , 'name' , 'epsilon_y' ;
};
M_.mapping.y.eqidx = [1 2 ];
M_.mapping.pi.eqidx = [1 2 ];
M_.mapping.epsilon_y.eqidx = [2 4 ];
M_.mapping.epsilon_pi.eqidx = [1 3 ];
M_.mapping.nu_y.eqidx = [4 ];
M_.mapping.nu_pi.eqidx = [3 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 2;
M_.block_structure.block(1).mfs = 2;
M_.block_structure.block(1).equation = [ 3 4];
M_.block_structure.block(1).variable = [ 4 3];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 4;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 2 3 4 ];
M_.block_structure.block(2).Simulation_Type = 8;
M_.block_structure.block(2).endo_nbr = 2;
M_.block_structure.block(2).mfs = 2;
M_.block_structure.block(2).equation = [ 2 1];
M_.block_structure.block(2).variable = [ 1 2];
M_.block_structure.block(2).is_linear = true;
M_.block_structure.block(2).NNZDerivatives = 7;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [1 3 4 6 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([1 1 2 1 2 1 2 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([1 3 3 4 4 6 6 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 2 2 4 6 6 8 ]);
M_.block_structure.variable_reordered = [ 4 3 1 2];
M_.block_structure.equation_reordered = [ 3 4 2 1];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 2 1;
 3 4;
 4 3;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 1;
 1 2;
 1 4;
 2 1;
 2 2;
 2 3;
 3 4;
 4 3;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 2;
 2 2;
];
M_.block_structure.dyn_tmp_nbr = 0;
M_.state_var = [4 3 1 ];
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
M_.NNZDerivatives = [15; 0; -1; ];
M_.dynamic_g1_sparse_rowval = int32([2 4 3 1 2 1 2 2 4 1 3 1 2 4 3 ]);
M_.dynamic_g1_sparse_colval = int32([1 3 4 5 5 6 6 7 7 8 8 10 10 13 14 ]);
M_.dynamic_g1_sparse_colptr = int32([1 2 2 3 4 6 8 10 12 12 14 14 14 15 16 ]);
M_.dynamic_g2_sparse_indices = int32([]);
M_.lhs = {
'pi'; 
'y'; 
'epsilon_pi'; 
'epsilon_y'; 
};
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 3];
M_.block_structure_stat.block(1).variable = [ 4];
M_.block_structure_stat.block(2).Simulation_Type = 3;
M_.block_structure_stat.block(2).endo_nbr = 1;
M_.block_structure_stat.block(2).mfs = 1;
M_.block_structure_stat.block(2).equation = [ 4];
M_.block_structure_stat.block(2).variable = [ 3];
M_.block_structure_stat.block(3).Simulation_Type = 6;
M_.block_structure_stat.block(3).endo_nbr = 2;
M_.block_structure_stat.block(3).mfs = 2;
M_.block_structure_stat.block(3).equation = [ 1 2];
M_.block_structure_stat.block(3).variable = [ 2 1];
M_.block_structure_stat.variable_reordered = [ 4 3 2 1];
M_.block_structure_stat.equation_reordered = [ 3 4 1 2];
M_.block_structure_stat.incidence.sparse_IM = [
 1 1;
 1 2;
 1 4;
 2 1;
 2 2;
 2 3;
 3 4;
 4 3;
];
M_.block_structure_stat.tmp_nbr = 0;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([1 2 1 2 ]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([1 1 2 2 ]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([1 3 5 ]);
M_.static_g1_sparse_rowval = int32([1 2 1 2 2 4 1 3 ]);
M_.static_g1_sparse_colval = int32([1 1 2 2 3 3 4 4 ]);
M_.static_g1_sparse_colptr = int32([1 3 5 7 9 ]);
load param_illus_det;
set_param_value('beta',beta); 
set_param_value('kappa',kappa); 
set_param_value('eta',eta); 
set_param_value('sigma',sigma); 
set_param_value('phi_pi',phi_pi); 
set_param_value('sigma_pi',sigma_pi); 
set_param_value('sigma_y',sigma_y); 
set_param_value('rho_pi',rho_pi); 
set_param_value('rho_y',rho_y); 
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(7))^2;
M_.Sigma_e(2, 2) = (M_.params(6))^2;
options_.drop = 2500;
options_.irf = 100;
options_.order = 1;
options_.periods = 200000;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'illus_det_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
