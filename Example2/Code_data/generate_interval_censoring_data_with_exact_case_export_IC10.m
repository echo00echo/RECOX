%% parameter setting 
%Target:
%IC: 10%
%Exact: 5%
%Right: 40%
%Left: 45%
IC = 10;  % number of different exam times 
rs = inf;  % truncation of right censor time
n = 500;   % samples per recovery
tau = 3;   % end time of sensing
% beta = [-1.2; zeros(48,1); -0.8; zeros(49,1); 0.6; zeros(49,1); 0.9; zeros(49,1); 1.5; zeros(800,1)]; % true beta
beta = [-1.2; zeros(48,1); -0.8; zeros(49,1); 0.6; zeros(49,1); 0.9; zeros(49,1); 1.5]; % true beta
p = length(beta);   % feature dimension

cov = 'Identity';	% covariance type
r = 0.5;   % covariance generation parameter

rep = 1; % replications

%% reset random seed
rng('default');

for rep_rr = 1:rep

%% generate mat file name


    generate_interval_censoring_data_with_exact_case(rep_rr, n, p, rs, tau, IC, cov, r, beta, rep);
end
    
function [Z, L, R, beta_tr] = generate_interval_censoring_data_with_exact_case(rep_rr, n, p, rs, tau, IC, cov, r, beta_tr, rep)
    mat_file_prefix_1 = sprintf('E_I%d_n%d_%s_rep%d', IC, n, cov, rep_rr);
    mat_filename = strcat(mat_file_prefix_1, '.mat');
    
% p = length(beta_tr);

switch cov
    case 'AR1'
        VarCovS = power(r,abs(repmat(1:p,p,1) - repmat((1:p).',1,p)));
    case 'CS'
        VarCovS = speye(p) + r*triu(ones(p,p),1) + r*tril(ones(p,p),-1);
    case 'Identity'
        VarCovS = speye(p);
    otherwise
        error('unknown input: cov');
end


Z = mvnrnd(zeros(p,1), VarCovS, n);

%% True failure times
% warning: Matlab's exprnd uses PDF f(x|mu) = 1/mu exp(-x/mu).  
T = 2*exprnd(exp(-Z*beta_tr));

%% number of exam times
CT = 10;

switch IC
    case 10
        timep = linspace(0,tau,CT);
        remove_ratio_interval = 0.4;      % ratio to remove from IC to LC
%       remove_ratio_right = 0.73;         % ratio to remove from LC to RC
        remove_ratio_exact = 0.33;         % ratio to remove from IC to exact
    case 30
        timep = linspace(0,tau,CT);
        remove_ratio = 1/2;         % ratio to remove
    case 50
        timep = linspace(0.00002,tau+20,CT);
    case 70  % newly added
        timep = linspace(0.00002,tau+7,CT);
    case 90
        timep = linspace(0.0008,tau+43,CT);
        remove_ratio = 1/3;         % ratio to remove
    otherwise
        error('unknown input: IC')
end

%% generate LR
if IC == 30 || IC == 50 || IC == 10
    obsT = repmat(timep, n, 1) .* binornd(1, 0.5, [n CT]);
elseif IC == 70  % newly added
    obsT = repmat(timep, n, 1) .* binornd(1, 0.8, [n CT]);
elseif IC == 90
    obsT = repmat(timep, n, 1) .* binornd(1, 0.95, [n CT]);
end

L = max(obsT.*(obsT < repmat(T,1,CT)), [], 2);

R = obsT.*(obsT > repmat(T,1,CT));
R(R==0) = inf;
R = min(R, [], 2);

if ~isempty(rs)
    R(R == inf) = rs;
end

%% assign censoring type
left_censoring = (L==0)&(T<tau);   % qualified left censoring indicator
right_censoring = (~left_censoring)&(R==rs);   % right censoring indicator
intvl_censoring = (~left_censoring)&(~right_censoring)&(T<tau); % qualified interval censoring indicator

%if IC ~= 90
%    make_exact = random_true_remover(left_censoring, remove_ratio);
%    L(make_exact) = T(make_exact);
%    R(make_exact) = T(make_exact);
%end

%make_exact = random_true_remover(intvl_censoring, remove_ratio);
%L(make_exact) = T(make_exact);
%R(make_exact) = T(make_exact);


if IC == 10
    [new_intvl_censoring, make_left] = random_true_remover(intvl_censoring, remove_ratio_interval, n);
    L(make_left) = 0;
    R(make_left) = T(make_left);  % move 40% from IC to LC 


    [not_used, make_exact] = random_true_remover(new_intvl_censoring, remove_ratio_exact, n);
    L(make_exact) = T(make_exact);
    R(make_exact) = T(make_exact);  % move 95% from IC to RC 
end


%% update censoring type
left_censoring_new = (L==0)&(T<tau);   % qualified left censoring indicator
right_censoring_new = (~left_censoring_new)&(R==rs);   % right censoring indicator
exact_case = (L==R);
intvl_censoring_new = (~exact_case)&(~left_censoring_new)&(~right_censoring_new)&(T<tau); % qualified interval censoring indicator

    %% save current results
%        save(mat_filename, ...
 %     'Z', 'L', 'R', 'beta_tr', 'cov', 'r', 'rep', 'tau', '-v7.3');
 
   save(mat_filename, ...
      'not_used', 'make_exact', 'make_left', 'new_intvl_censoring','Z', 'L', 'R', 'beta_tr', 'cov', 'r', 'rep', 'tau', 'left_censoring', 'right_censoring', 'intvl_censoring', 'exact_case','left_censoring_new', 'right_censoring_new', 'intvl_censoring_new');
  
% 
% %% original 
% for i = 1: length(R)
%     if L(i) == 0
%         ctype(i)=1;
%     elseif R(i) == rs
%         ctype(i)=2;
%     else
%         ctype(i)=3;
%     end
% end
% 
% %% assign exact case
% data=sortrows([L R T ctype.'], 4);
% 
% %% from left censoring
% p = randperm(length(ctype(ctype==1)));
% p1 = p(1:floor(length(p)/2));
% 
% % set new values.
% data(p1,1)=data(p1, 3);
% data(p1,2)=data(p1, 3);
% 
% %% from interval censoring
% p3 = randperm(length(ctype(ctype==3)));
% p4 = p3(1:floor(length(p3)/2)) + length((ctype(ctype==1))) + length((ctype(ctype==2)));
% 
% % set new values.
% data(p4,1)=data(p4, 3);
% data(p4,2)=data(p4, 3);
% 
% %% update censorying type
% for i = 1: length(R)
%     if data(i, 1) == 0
%         data(i, 4)=1;
%     elseif data(i, 2) == rs
%         data(i, 4)=2;
%     elseif data(i, 1)==data(i, 2)
%         data(i, 4)=3;
%     else
%         data(i, 4)=4;
%     end
% end
% 
% L = data(:, 1);
% R = data(:, 2);
% ctype = data(:, 4);
% % LR = [L R];
% % if ~isempty(rs)
% %     LR(LR == inf) = rs;
% % end

end

%% given a logical array, return a logical array with a randomly selected subset of true's
function [ind_vec, remove] = random_true_remover(ind_vec, remove_ratio, n)
    remove = false(n);
    loc = find(ind_vec);    % position indices of true's in ind_vec
    ntrue = numel(loc);     % the number of true's
    idx_to_remov = randperm(ntrue, floor(ntrue*remove_ratio));    % randomly select some from 1:ntrue
    ind_vec(loc(idx_to_remov)) = false;
    remove(loc(idx_to_remov)) = true;
end
