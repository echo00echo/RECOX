function generate_failure_times_and_interval_censoring_data
% This function generates synthetic data for cox-model regression

%% parameter setting
rs = inf;   % truncation of right censor time, either a finite value or inf (default is inf)
n = 1000;   % number of samples
T_end = 3;  % study end time; study always starts at time 0
n_exams = 10;  % number of exams between time 0 and T_end
showup_rate = 0.5; % the probability each individual showing up at each exam
IC = 10; % interval censoring ratio

% true beta, coefficients of the feature vector
beta_tr = [-1; zeros(48,1); -0.8; zeros(49,1); 0.6; zeros(49,1); 0.9; zeros(49,1); 1.5; zeros(4300,1)]; 
p = length(beta_tr);  % feature dimension

cov = 'AR1';    % covariance structure of features: 'AR1' = auto-regressive, order 1; 'CS' = compound symmetric; 'Identity' 
r = 0.5;        % compound symmetric's parameter, only for cov = 'CS'

rep = 100; % replications

% reset random seed
rng('default');



%% generate data and save to mat-files
for rep_rr = 1:rep
    generate_interval_censoring_data_with_exact_case(rep_rr, IC, n, p, rs, n_exams, showup_rate, T_end, cov, r, beta_tr, rep);
    

 % generate matfile names
mat_file_prefix = sprintf('E_I%d_n%d_%s_rep%d', IC, n, cov, rep_rr);
mat_filename = strcat(mat_file_prefix, '.mat');

% data loading preparation and load cenroing ratios for each type
data_folder = pwd;

    left_c = load(strcat(data_folder, '\', mat_filename), 'left_count');
    left_r(rep_rr) = vertcat(left_c.left_count)*100/n;

    right_c = load(strcat(data_folder, '\', mat_filename), 'right_count');
    right_r(rep_rr) = vertcat(right_c.right_count)*100/n;

    intvl_c = load(strcat(data_folder, '\', mat_filename), 'intvl_count');
    intvl_r(rep_rr) = vertcat(intvl_c.intvl_count)*100/n;

    exact_c = load(strcat(data_folder, '\', mat_filename), 'exact_count');
    exact_r(rep_rr) = vertcat(exact_c.exact_count)*100/n;

end
    save(strcat(mat_file_prefix, '_ratio', '.mat'), 'left_r', 'right_r', 'intvl_r', 'exact_r');
end % end of main function


function [Z, L, R, beta_tr] = generate_interval_censoring_data_with_exact_case(rep_idx, IC, n, p, rs, n_exams, showup_rate, T_end, cov, r, beta_tr, reps)

% truncation of right censor time
if isempty(rs); rs = inf; end

% generate matfile names
mat_file_prefix = sprintf('E_I%d_n%d_%s_rep%d', IC, n, cov, rep_idx);
mat_filename = strcat(mat_file_prefix, '.mat');
    
%% generate feature matrix Z
switch cov
    case 'AR1'  % auto-regressive, order 1
        VarCovS = power(r,abs(repmat(1:p,p,1) - repmat((1:p).',1,p)));
    case 'CS'   % compound symmetric
        VarCovS = speye(p) + r*triu(ones(p,p),1) + r*tril(ones(p,p),-1);
    case 'Identity'
        VarCovS = speye(p);
    otherwise
        error('unknown input: cov');
end

Z = mvnrnd(zeros(p,1), VarCovS, n);


%% generate n i.i.d. samples of L and R from true failure time
% [L, R] = generate_L_R(n, Z, beta_tr, rs, n_exams, showup_rate, T_end);

% Sample true failure times
T = exprnd(exp(-Z*beta_tr));
% Note: this generation corresponds to a constant hazard rate for each individual
% for individual i, lam_i(t) = exp(-Z_i'*beta_tr)*lam0(t) with baseline rate lam0(t) = 1

% show up at any exam at the rate of showup_rate
shows = binornd(1, showup_rate, [n n_exams]);

% re-sample those who never show up
never_show = (sum(shows,2)==0); % the number of showups of individuals
while any(never_show)
    shows(never_show,:) = binornd(1, showup_rate, [nnz(never_show) n_exams]);
    never_show = (sum(shows,2)==0);
end

% equally spaced exam times in (0, T_end], where t=0 is not an exam time.
exam_times = linspace(0, T_end, n_exams+1);
exam_times = exam_times(2:end);

% observation time matrix, each 0-entry is a no-show
obs_times = repmat(exam_times, n, 1) .* shows;

% compute the last time when each individual shows up alive 
L = max(obs_times.*(obs_times < repmat(T,1,n_exams)), [], 2);  % 2 means taking max along dim 2

% compute the first time when each individual "shows up" dead
R = obs_times.*(obs_times > repmat(T,1,n_exams));
R(R==0) = inf;  % never found dead
R = min(R, [], 2);

if ~isempty(rs) && ~isinf(rs)
    R(R == inf) = rs;
end

%% when the generated left, right, intvl, and exact censoring cases are not ideal, redistribute the cases

assert(all((L>0)|(R<rs)), 'found sample(s) with L==0 and R==rs, impossible');

% generate logical arrays
left_censoring = (L==0);                        % left censoring indicator
right_censoring = (~left_censoring)&(R==rs);    % right censoring indicator
intvl_censoring = (~left_censoring)&(~right_censoring)&(L<R); % interval censoring indicator
exact_censoring = (~left_censoring)&(~right_censoring)&(~intvl_censoring)&(L==R); % exact censoring indicator

fprintf('\nrep %2d       left      |    right     |    intvl     |    exact    \n', rep_idx);
fprintf('original: %4d (%4.1f%%) | %4d (%4.1f%%) | %4d (%4.1f%%) | %4d (%4.1f%%) \n', ...
    nnz(left_censoring), nnz(left_censoring)*100/n, ...
    nnz(right_censoring), nnz(right_censoring)*100/n, ...
    nnz(intvl_censoring), nnz(intvl_censoring)*100/n, ...
    nnz(exact_censoring), nnz(exact_censoring)*100/n );

% left censoring --> exact censoring
move_left_to_exact = 0.0; 
if move_left_to_exact > 0
    [left_censoring, make_exact] = random_true_remover(left_censoring, move_left_to_exact);
    L(make_exact) = min(T(make_exact), T_end);
    R(make_exact) = min(T(make_exact), T_end);
end

% left censoring --> right censoring
move_left_to_right = 1;
if move_left_to_right > 0
    [left_censoring, make_right] = random_true_remover(left_censoring, move_left_to_right);
    L(make_right) = min(T(make_right), T_end);
    R(make_right) = inf;
end

% intvl censoring --> right censoring
move_intvl_to_right = 0.6;
if move_intvl_to_right > 0
    [intvl_censoring, make_right] = random_true_remover(intvl_censoring, move_intvl_to_right);
    L(make_right) = T(make_right);
    R(make_right) = inf;
end



%% update censoring type
assert(isequal(left_censoring, L==0));   % qualified left censoring indicator
right_censoring = (~left_censoring)&(R==rs);   % right censoring indicator
exact_censoring = (L==R);
assert(isequal(intvl_censoring, (~exact_censoring)&(~left_censoring)&(~right_censoring))); % interval censoring indicator

left_count = nnz(left_censoring);
right_count = nnz(right_censoring);
intvl_count = nnz(intvl_censoring); 
exact_count = nnz(exact_censoring); 

fprintf('adjusted: %4d (%4.1f%%) | %4d (%4.1f%%) | %4d (%4.1f%%) | %4d (%4.1f%%) \n\n', ...
    nnz(left_censoring), nnz(left_censoring)*100/n, ...
    nnz(right_censoring), nnz(right_censoring)*100/n, ...
    nnz(intvl_censoring), nnz(intvl_censoring)*100/n, ...
    nnz(exact_censoring), nnz(exact_censoring)*100/n );
 
save(mat_filename, 'n', 'Z', 'T', 'L', 'R', 'beta_tr', 'cov', 'r', 'rep_idx', 'reps', 'T_end', 'left_censoring', 'right_censoring', 'intvl_censoring', 'exact_censoring', 'left_count', "right_count", "intvl_count", "exact_count");
  
end


%% Given a logical array 'ind_vec', return 'ind_vec' with approx 'remove_ratio' 
% of true's in 'ind_vec', randomly selected, changed to false's.
% 'flip_loc' indicates the components of ind_vec that have been changed 

function [ind_vec, flip_loc] = random_true_remover(ind_vec, remove_ratio)

    true_loc = find(ind_vec);    % position indices of true's in ind_vec
    ntrue = numel(true_loc);     % the number of true's

    idx_to_remov = randperm(ntrue, floor(ntrue*remove_ratio));    % randomly select a subset from 1:ntrue
    ind_vec(true_loc(idx_to_remov)) = false;

    flip_loc = false(size(ind_vec));
    flip_loc(true_loc(idx_to_remov)) = true;
end
