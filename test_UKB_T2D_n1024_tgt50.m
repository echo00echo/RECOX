%% parameter setting
IC = 1; 
n = 1024;   % samples per recovery
p = 5201;   % feature dimension
cov = 'UKB_T2D';	% covariance type
% r = [];   % covariance generation parameter

sp = 1; % starting point of # replications in case used in more than 1 server
rep = 1; % replications

B = 1000;  % rounds of random splitting 

%% call test function
generic_opts.estimation_warmstart = true;
generic_opts.use_parfor = true;
generic_opts.target_nnz = 50;
generic_opts.real_data = true;
generic_opts.chr = 3910; 
% opts.RandSeed = 2021;

generic_opts.mult_suff_j = true;
generic_opts.mult_suff_j_blk = 5;
generic_opts.rngseed = 11;

screen_opts.regularization = 'l1';  % 'none', 'l2 square', 'l1'
screen_opts.BCDBlockSize = 10;
screen_opts.BCDEpochNum = 10;
screen_opts.GreedyInterval = 20;
screen_opts.GreedyRule = 'Gauss-Southwell-r ls';
screen_opts.DisplayLevel = 5; % 1

restr_est_opts.regularization = 'none'; % 'none', 'l2 square', 'l1'
restr_est_opts.theta = 1e-3;
restr_est_opts.kn_beta_bnd = 10; 
restr_est_opts.DisplayLevel = 1;

% test_splitting_and_smoothing_multi_rep_3(IC, n, p, cov, rep, B, generic_opts, screen_opts, restr_est_opts);
test_cvx_splitting_and_smoothing_2(IC, n, p, cov, sp, rep, B, generic_opts, screen_opts, restr_est_opts);

