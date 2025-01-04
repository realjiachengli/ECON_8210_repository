%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clearvars -global
clear_persistent_variables(fileparts(which('dynare')), false)
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'perturb_dynare';
M_.dynare_version = '6.2';
oo_.dynare_version = '6.2';
options_.dynare_version = '6.2';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'e'};
M_.exo_names_tex(1) = {'e'};
M_.exo_names_long(1) = {'e'};
M_.endo_names = cell(4,1);
M_.endo_names_tex = cell(4,1);
M_.endo_names_long = cell(4,1);
M_.endo_names(1) = {'c'};
M_.endo_names_tex(1) = {'c'};
M_.endo_names_long(1) = {'c'};
M_.endo_names(2) = {'l'};
M_.endo_names_tex(2) = {'l'};
M_.endo_names_long(2) = {'l'};
M_.endo_names(3) = {'k'};
M_.endo_names_tex(3) = {'k'};
M_.endo_names_long(3) = {'k'};
M_.endo_names(4) = {'z'};
M_.endo_names_tex(4) = {'z'};
M_.endo_names_long(4) = {'z'};
M_.endo_partitions = struct();
M_.param_names = {};
M_.param_names_tex = {};
M_.param_names_long = {};
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 4;
M_.param_nbr = 0;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
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
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.nonzero_hessian_eqs = [2 3 4];
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
 0 3 7;
 0 4 8;
 1 5 0;
 2 6 9;]';
M_.nstatic = 0;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [9; 11; 16; 4; ];
M_.equations_tags = {
  1 , 'name' , 'z' ;
  2 , 'name' , 'lam' ;
  3 , 'name' , 'l' ;
  4 , 'name' , 'k' ;
};
M_.mapping.c.eqidx = [2 3 4 ];
M_.mapping.l.eqidx = [2 3 4 ];
M_.mapping.k.eqidx = [2 3 4 ];
M_.mapping.z.eqidx = [1 2 3 4 ];
M_.mapping.e.eqidx = [1 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 1;
M_.block_structure.block(1).mfs = 1;
M_.block_structure.block(1).equation = [ 1];
M_.block_structure.block(1).variable = [ 4];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 2;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 2 ];
M_.block_structure.block(2).Simulation_Type = 8;
M_.block_structure.block(2).endo_nbr = 3;
M_.block_structure.block(2).mfs = 3;
M_.block_structure.block(2).equation = [ 4 2 3];
M_.block_structure.block(2).variable = [ 3 1 2];
M_.block_structure.block(2).is_linear = false;
M_.block_structure.block(2).NNZDerivatives = 11;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [1 4 5 6 8 9 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([1 1 2 3 1 2 3 1 3 2 2 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([1 4 4 4 5 5 5 6 6 8 9 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 2 2 2 5 8 10 10 11 12 ]);
M_.block_structure.variable_reordered = [ 4 3 1 2];
M_.block_structure.equation_reordered = [ 1 4 2 3];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 1 4;
 4 3;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 4;
 2 1;
 2 3;
 3 1;
 3 2;
 3 3;
 3 4;
 4 1;
 4 2;
 4 3;
 4 4;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 2 1;
 2 2;
 2 4;
];
M_.block_structure.dyn_tmp_nbr = 8;
M_.state_var = [4 3 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(0, 1);
M_.endo_trends = struct('deflator', cell(4, 1), 'log_deflator', cell(4, 1), 'growth_factor', cell(4, 1), 'log_growth_factor', cell(4, 1));
M_.NNZDerivatives = [17; 42; 51; ];
M_.dynamic_g1_sparse_rowval = int32([4 1 2 3 4 3 4 2 3 4 1 3 4 2 2 2 1 ]);
M_.dynamic_g1_sparse_colval = int32([3 4 5 5 5 6 6 7 7 7 8 8 8 9 10 12 13 ]);
M_.dynamic_g1_sparse_colptr = int32([1 1 1 2 3 6 8 11 14 15 16 16 17 18 ]);
M_.dynamic_g2_sparse_indices = int32([2 5 5 ;
2 9 9 ;
2 9 10 ;
2 9 7 ;
2 9 12 ;
2 10 10 ;
2 10 7 ;
2 10 12 ;
2 7 7 ;
2 7 12 ;
2 12 12 ;
3 5 5 ;
3 5 6 ;
3 5 7 ;
3 5 8 ;
3 6 6 ;
3 6 7 ;
3 6 8 ;
3 7 7 ;
3 7 8 ;
3 8 8 ;
4 6 6 ;
4 6 3 ;
4 6 8 ;
4 3 3 ;
4 3 8 ;
4 8 8 ;
]);
M_.dynamic_g3_sparse_indices = int32([2 5 5 5 ;
2 9 9 9 ;
2 9 9 10 ;
2 9 9 7 ;
2 9 9 12 ;
2 9 10 10 ;
2 9 10 7 ;
2 9 10 12 ;
2 9 7 7 ;
2 9 7 12 ;
2 9 12 12 ;
2 10 10 10 ;
2 10 10 7 ;
2 10 10 12 ;
2 10 7 7 ;
2 10 7 12 ;
2 10 12 12 ;
2 7 7 7 ;
2 7 7 12 ;
2 7 12 12 ;
2 12 12 12 ;
3 5 5 5 ;
3 5 5 6 ;
3 5 5 7 ;
3 5 5 8 ;
3 5 6 6 ;
3 5 6 7 ;
3 5 6 8 ;
3 5 7 7 ;
3 5 7 8 ;
3 5 8 8 ;
3 6 6 6 ;
3 6 6 7 ;
3 6 6 8 ;
3 6 7 7 ;
3 6 7 8 ;
3 6 8 8 ;
3 7 7 7 ;
3 7 7 8 ;
3 7 8 8 ;
3 8 8 8 ;
4 6 6 6 ;
4 6 6 3 ;
4 6 6 8 ;
4 6 3 3 ;
4 6 3 8 ;
4 6 8 8 ;
4 3 3 3 ;
4 3 3 8 ;
4 3 8 8 ;
4 8 8 8 ;
]);
M_.lhs = {
'z'; 
'1/c'; 
'l'; 
'k'; 
};
M_.static_tmp_nbr = [8; 2; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 1];
M_.block_structure_stat.block(1).variable = [ 4];
M_.block_structure_stat.block(2).Simulation_Type = 6;
M_.block_structure_stat.block(2).endo_nbr = 3;
M_.block_structure_stat.block(2).mfs = 3;
M_.block_structure_stat.block(2).equation = [ 2 3 4];
M_.block_structure_stat.block(2).variable = [ 3 2 1];
M_.block_structure_stat.variable_reordered = [ 4 3 2 1];
M_.block_structure_stat.equation_reordered = [ 1 2 3 4];
M_.block_structure_stat.incidence.sparse_IM = [
 1 4;
 2 1;
 2 2;
 2 3;
 2 4;
 3 1;
 3 2;
 3 3;
 3 4;
 4 1;
 4 2;
 4 3;
 4 4;
];
M_.block_structure_stat.tmp_nbr = 7;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 2 3 1 2 3 1 2 3 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 1 1 2 2 2 3 3 3 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 4 7 10 ]);
M_.static_g1_sparse_rowval = int32([2 3 4 2 3 4 2 3 4 1 2 3 4 ]);
M_.static_g1_sparse_colval = int32([1 1 1 2 2 2 3 3 3 4 4 4 4 ]);
M_.static_g1_sparse_colptr = int32([1 4 7 10 14 ]);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(4) = 0;
oo_.steady_state(3) = 3.7612;
oo_.steady_state(2) = 0.9465;
oo_.steady_state(1) = 1.1161;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
steady;
options_.k_order_solver = true;
options_.order = 3;
options_.pruning = true;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'perturb_dynare_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
