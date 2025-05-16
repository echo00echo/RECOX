%% parameter setting
IC = 10; 
n = 1000;   % samples per recovery
p = 4500;   % feature dimension
cov = 'CS';	% covariance type
r = 0.05;   % covariance generation parameter

rep = 100; % replications

B = 1000;  % rounds of random splitting 

%% call test function
generic_opts.estimation_warmstart = true;
generic_opts.use_parfor = true;
generic_opts.target_nnz = n/log(n);
generic_opts.gen_exact = true;
%generic_opts.chr = 6; 
% opts.RandSeed = 2021;
%generic_opts.resume = true;
%generic_opts.resume_at_rr = 3;


generic_opts.mult_suff_j = true;
generic_opts.mult_suff_j_blk = 5;
generic_opts.rngseed = 1100;

screen_opts.regularization = 'l1';  % 'none', 'l2 square', 'l1'
screen_opts.BCDBlockSize = 10;
% screen_opts.BCDEpochNum = 10;
% screen_opts.GreedyInterval = 20;
screen_opts.GreedyRule = 'Gauss-Southwell-r ls';
screen_opts.DisplayLevel = 1; % 2




restr_est_opts.regularization = 'none'; % 'none', 'l2 square', 'l1'
restr_est_opts.theta = 1e-3;
restr_est_opts.beta_bnd = 10; 
restr_est_opts.DisplayLevel = 1;

% test_splitting_and_smoothing_multi_rep_3(IC, n, p, cov, rep, B, generic_opts, screen_opts, restr_est_opts);
test_cvx_splitting_and_smoothing_read_p_resume(IC, n, p, cov, rep, B, generic_opts, screen_opts, restr_est_opts);
