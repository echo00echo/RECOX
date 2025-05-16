function test_cvx_splitting_and_smoothing_read_p(IC, n, p, cov, rep, B, generic_opts, screen_opts, restr_est_opts)

clc
close all;
% clear;

%% Set default options and then read the input options
opts.chr = [];
opts.RandSeed = [];
opts.use_parfor = true;
opts.estimation_warmstart = true;
opts.target_nnz = min(round(n/log(n)),p);
opts.gen_exact = false;
opts.real_data = false;
opts.reduced_data = false;
% ------------ for RegParameter search ------------- 
opts.regpar.Low = 5e-2; % lower bound for RegParameter search
opts.regpar.Hgh = 5e1;  % upper bound for RegParameter search 
opts.regpar.Tol = 1e-1; % stop tolerance for RegParameter search
% ------------ for multiple_shuffled_j in partial regression ------------- 
opts.mult_suff_j = false;
opts.mult_suff_j_blk = 5; % number of j indices in each partial regression
% ------------ seed for Randperm ------------- 
opts.rngseed = [];
% ------------ for resume the for-loop ------------- 
opts.resume = false;
  opts.resume_at_rr = [];
  opts.filename = [];

if nargin >= 7 && ~isempty(generic_opts)
    opts = ReadOptions(opts, generic_opts);
    if opts.target_nnz > p
        opts.target_nnz = p;
    end
end

%% generate mat file name
mat_file_prefix = sprintf('I%d_n%d_p%d_%s_tgt%d_rep%d_B%d_block%d',IC,n,p,cov,opts.target_nnz,rep,B, opts.mult_suff_j_blk);

%% data loading preparation and load beta_tr
data_folder = strcat(fileparts(which(mfilename)), '/../originaldata/');
%data_folder = strcat(fileparts(which(mfilename)), 'C:/Users/echo0/Dropbox/Share_Rui/IC/repo/originaldata/');

rr = 1;
if opts.gen_exact
    data_filename = sprintf('E_I%d_n%d_%s_rep%d', IC, n, cov, rr);
elseif opts.real_data
        data_filename = sprintf('E_IC%d_n%d_%s_chr%d_p%d', IC, n, cov, opts.chr, p);
elseif opts.reduced_data
        data_filename = sprintf('R_IC%d_n%d_%s_chr%d', IC, n, cov, opts.chr);
else
    data_filename = sprintf('I%d_n%d_%s_rep%d', IC, n, cov, rr);
end
load(strcat(data_folder, data_filename, '.mat'), 'beta_tr');


%% get true support and nnz
if exist('beta_tr','var') && ~isempty(beta_tr) 
    beta_tr = beta_tr(1:p);
    S_tr = find(beta_tr);
    nz = nnz(beta_tr);
    if nz == 0; error('True beta is not given in beta_tr'); end
else
    beta_tr = [];
end

%% initialize variables to save
rnd_strm_start = RandStream.getGlobalStream;
beta_supp_est_save = zeros(p,B,rep);
beta_save = zeros(p,B,rep);
loss_save = zeros(p,B,rep);
samp_save = false(n,B,rep);
if ~isempty(beta_tr) 
    betO_save = zeros(nz,rep);
end
rand_save = zeros(numel(get(rnd_strm_start,'state')),rep);

%% options for beta's support estimation
supp_est_opts.DisplayLevel = 1;
supp_est_opts.Method = 'alternating minimization by greedy BCGD';
    % 'joint minimization by solver'
    % 'alternating minimization by solver'
    % 'alternating minimization by BCGD'
    % 'alternating minimization by greedy BCGD'
supp_est_opts.OuterItr = 40;
supp_est_opts.regularization = 'l1';  % 'none', 'l2 square', 'l1'
supp_est_opts.BCDBlockSize = 10;
supp_est_opts.BCDEpochNum = 4;
supp_est_opts.GreedyInterval = 5;
supp_est_opts.GreedyRule = 'Gauss-Southwell-r ls';
    % 'noreg gradient'
    % 'Gauss-Southwell-r'
    % 'Gauss-Southwell-s'
supp_est_opts.DispInnerInterval = 100;
%est_opts.InitialBeta = [];
if isfield(opts, 'RandSeed')
    supp_est_opts.RandSeed = opts.RandSeed;
end
supp_est_opts.TrueSolution = beta_tr; % used for progress print

if nargin >= 8 && ~isempty(screen_opts)
    supp_est_opts = ReadOptions(supp_est_opts, screen_opts);
end


%% options for restricted estimation
res_est_opts.DisplayLevel = 1;
res_est_opts.DispInnerInterval = 500;
res_est_opts.regularization = 'none';
res_est_opts.theta = 1e-3;
res_est_opts.beta_bnd = 2;

if nargin >= 9 && ~isempty(restr_est_opts)
    res_est_opts = ReadOptions(res_est_opts, restr_est_opts);
end


%% options for oracle estimate
orcl_opts.DisplayLevel = 0.5;
orcl_opts.TrueSolution = beta_tr;

%% start diary
if opts.resume && ~isempty(opts.filename) 
    diary(strcat(opts.filename,'.log'));
    mat_filename = strcat(opts.filename,'.mat');
else
    dattim_str = datestr(datetime('now'),'yyyymmmddTHHMM');
    diary(strcat(mat_file_prefix,'_',dattim_str,'.log'));
    mat_filename = strcat(mat_file_prefix,'_',dattim_str,'.mat');
end
    
%% main replication loop
for rr = 2:rep
    
    % only used for resuming from rr-th repetition
    if opts.resume && rr < opts.resume_at_rr - 1
        continue;
    elseif opts.resume && rr == opts.resume_at_rr - 1
        fprintf('load %s at rr = %d\n\n', mat_filename, rr);
        load temp.mat;
        % set(rnd_strm_start, 'state', rand_save(:,rr));
    else
        % NORMAL mode
        % save the random generator's state at the beginning of a repl loop
        rand_save(:,rr) = get(rnd_strm_start, 'state');
    end
    
    %% data loading
    if opts.gen_exact
        data_filename = sprintf('E_I%d_n%d_%s_rep%d', IC, n, cov, rr);
    elseif opts.real_data
        data_filename = sprintf('E_IC%d_n%d_%s_chr%d_p%d', IC, n, cov, opts.chr, p);
    elseif opts.reduced_data
        data_filename = sprintf('R_IC%d_n%d_%s_chr%d', IC, n, cov, opts.chr);
    else
        data_filename = sprintf('I%d_n%d_%s_rep%d', IC, n, cov, rr);
    end
    load(strcat(data_folder, data_filename, '.mat'), 'Z', 'L', 'R');
    Z = Z(:, 1:p);

    %% correct R, which is actually unnecessary as we support inf
%     if max(L) > 10 && max(L) < 30
%         R(R==rs) = 30;
%     elseif max(L) >= 30 && max(L) < 50
%         R(R==rs) = 50;
%     elseif max(L) >= 50
%        error('max(L) too large'); 
%     end

    %% compute oracle beta estimate
    if ~isempty(beta_tr)
        loss = intvl_cnsr_cvx_loglikelihood_class_1(size(Z,1), L, R);
        optz = intvl_cnsr_cvx_loglikelihood_maximize_class_1(Z(:,S_tr), loss, orcl_opts);
        try
            [~, ~, betO] = maximize_by_knitro(optz, []);
            betO_save(:,rr) = betO;
        catch ME
            clear mex;
            betO_save(:,rr) = inf;
            fprintf('\n\nError message: %s. rr=%d\n', ME.message, rr);
            keyboard;
        end
    end
    
    %% Main iteration of splitting
    % hbeta = zeros(p,1);
    for bb = 1:B
        % fprintf('Splitting round %d/%d:\n', bb, B); drawnow;
        
        %% random split
        rng(opts.rngseed + bb); % reset seed
        sample_perm = randperm(n);
        D_1 = sample_perm(1:floor(n/2));
        D_2 = sample_perm(floor(n/2)+1:end);
        Z_2 = Z(D_2,:);
        L_2 = L(D_2); R_2 = R(D_2);
        
        %% call solver to obtain S
        supp_est_opts.print.head = '  b type';
        supp_est_opts.print.text = sprintf('%3d scrn', bb);
        % try
            [phi_est, beta_supp_est] = support_estimation(L_2, R_2, Z_2, opts.target_nnz, supp_est_opts, opts.regpar);
            S = find(beta_supp_est);
            if ~isempty(beta_tr)
                fprintf('\nb=%d: miss=%d false=%d |beta-true|/|true|=%4.2e\n\n',...
                    bb, length(setdiff(S_tr, S)), length(setdiff(S, S_tr)),...
                    norm(beta_tr - beta_supp_est)/norm(beta_tr));
            else
                fprintf('\nb=%d: support estimation complete\n\n', bb);                
            end
%         catch ME
%             clear mex;
%             fprintf('\n\nsupport_estimation error: %s\nrr=%d, bb=%d\n', ME.message, rr, bb);
%             keyboard;
%         end


        %% compute beta tilde
    
        res_est_opts.print.head = {'  b type   j    i     obj     phi      beta '};
        fprintf('\n%s\n%s\n', res_est_opts.print.head{1}, repmat('-',1,length(res_est_opts.print.head{1})));
    
        res_est_opts.print.text = sprintf('% 3d estm', bb);
        % existing_len = length(res_opts.print.text);
    
        Z_1 = Z(D_1,:);
        L_1 = L(D_1); R_1 = R(D_1);
       
        %% save data splitting 
         split_filename = sprintf('E_I%d_n%d_%s_rep%d_split%d', IC, n, cov, rr, bb);
         save(split_filename, ...
        'IC', 'n', 'p', 'cov', ...  % problem parameters
        'rr', 'rep', 'bb', 'B', ... % current progress indicators
        'L_1', 'L_2', 'R_1', 'R_2', 'D_1', 'D_2');   % current split
        % Initialize tbeta vector
        tbeta = zeros(p,1);
        tloss = zeros(p,1);   % vector used to save loss values
        
        if opts.estimation_warmstart
            beta_S_est = beta_supp_est(S);
            if ~isempty(setxor([L_1 R_1], [L_2 R_2]))
                phi_est = [];   % avoid using the screening phi to warm start
            end
        else
            beta_S_est = [];
            phi_est = [];
        end
        
        % replace: tbeta(S) = beta_supp_est(S); 
        [temp, loss, phi_est] = restricted_estimation(L_1, R_1, Z_1, S, res_est_opts, [], beta_S_est, phi_est);
        tbeta(S) = temp;
        tloss(S) = loss;
        
        % save/reset beta_S_est and phi_est for warmstarting in restricted_estimation()
        if opts.estimation_warmstart
            beta_S_est = temp;
        else
            beta_S_est = [];
            phi_est = [];
        end
        
        j_set = setdiff(1:p, S);    % set of indices to run over estimation
        
        % Generate j_cells, cell array of index blocks
        if opts.mult_suff_j
            j_len = numel(j_set);
            j_set = j_set(randperm(j_len));
            j_cell_len = ceil(j_len/opts.mult_suff_j_blk);
            j_cells = cell(1,j_cell_len);
            for jj = 1:j_cell_len
                j_cells{jj} = sort(j_set((jj-1)*opts.mult_suff_j_blk+1 : min(jj*opts.mult_suff_j_blk,j_len)));
            end
        else
            j_cells = num2cell(j_set);
        end
        j_cell_len = length(j_cells);
        clear j_set j_len

        % Prepare for the for loop
        S_len = numel(S);
        if opts.mult_suff_j
            if opts.use_parfor
                % due to parfor restriction, set up extra variables
                parfor_tbeta = cell(1,j_cell_len);
                parfor_tloss = zeros(j_cell_len, 1);
                % tic;
                parfor jj = 1:j_cell_len
                    [temp, val] = restricted_estimation(L_1, R_1, Z_1, union(S,j_cells{jj},'stable'), res_est_opts, j_cells{jj}, beta_S_est, phi_est);  % TODO j
                    parfor_tbeta{jj} = temp(S_len+1:end);
                    parfor_tloss(jj) = val;
                end
                % toc;

                for jj = 1:j_cell_len
                    tbeta(j_cells{jj}) = parfor_tbeta{jj};
                    tloss(j_cells{jj}) = parfor_tloss(jj);
                end
            else
                for jj = 1:j_cell_len
                    [temp, val] = restricted_estimation(L_1, R_1, Z_1, union(S,j_cells{jj},'stable'), res_est_opts, j_cells{jj}, beta_S_est, phi_est);  % TODO j
                    tbeta(j_cells{jj}) = temp(S_len+1:end);
                    tloss(j_cells{jj}) = val;
                end
            end
        else
            if opts.use_parfor
                % due to the parfor restriction, set up extra variables
                parfor_tbeta = zeros(j_cell_len,1);
                parfor_tloss = zeros(j_cell_len,1);
                % tic;
                parfor jj = 1:j_cell_len
                    [temp, val] = restricted_estimation(L_1, R_1, Z_1, union(S,j_cells{jj},'stable'), res_est_opts, j_cells{jj}, beta_S_est, phi_est);  % TODO j
                    parfor_tbeta(jj) = temp(end);
                    parfor_tloss(jj) = val;
                end
                % toc;

                for jj = 1:j_cell_len
                    tbeta(j_cells{jj}) = parfor_tbeta(jj);
                    tloss(j_cells{jj}) = parfor_tloss(jj);
                end
            else
                for jj = 1:j_cell_len
                    [temp, val] = restricted_estimation(L_1, R_1, Z_1, union(S,j_cells{jj},'stable'), res_est_opts, j_cells{jj}, beta_S_est, phi_est);  % TODO j
                    tbeta(j_cells{jj}) = temp(end);
                    tloss(j_cells{jj}) = val;
                end
            end
        end
    %     fprintf('b=%d: miss=%d false=%d |beta-true|=%4.2e\n\n',...
    %             bb, length(setdiff(S_tr, find(tbeta))), length(setdiff(find(tbeta),S_tr)),...
    %             norm(beta_tr - tbeta)/norm(beta_tr));            
    %     hbeta = hbeta + tbeta;

        %% update storage
        beta_supp_est_save(:,bb,rr) = beta_supp_est;
        beta_save(:,bb,rr) = tbeta;
        loss_save(:,bb,rr) = tloss;
        samp_save(D_1,bb,rr) = true;

        %% save current results
        data_filename_b = sprintf('I%d_n%d_%s_rep%d_B%d_block%d', IC, n, cov, rr, bb, opts.mult_suff_j_blk);
    if ~isempty(beta_tr)
         save(data_filename_b, ...
        'IC', 'n', 'p', 'cov', ...  % problem parameters
        'rr', 'rep', 'bb', 'B', ... % current progress indicators
        'rand_save', ... % 'rnd_strm_start', 'rnd_state_last', ... % current rand gen state
        'beta_supp_est_save', 'beta_save', 'loss_save', 'samp_save', 'betO_save', ...   % current results
        '-v7.3'); % save data larger than 2 GB use MAT-file version 7.3 or later

    else
         save(data_filename_b, ...
        'IC', 'n', 'p', 'cov', ...  % problem parameters
        'rr', 'rep', 'bb', 'B', ... % current progress indicators
        'rand_save', ... % 'rnd_strm_start', 'rnd_state_last', ... % current rand gen state
        'beta_supp_est_save', 'beta_save', 'loss_save', 'samp_save', ...   % current results
        '-v7.3'); % save data larger than 2 GB use MAT-file version 7.3 or later
    end

    end
    
    %% save current results
    if ~isempty(beta_tr)
    save(mat_filename, ...
        'IC', 'n', 'p', 'cov', ...  % problem parameters
        'rr', 'rep', 'bb', 'B', ... % current progress indicators
        'rand_save', ... % 'rnd_strm_start', 'rnd_state_last', ... % current rand gen state
        'beta_supp_est_save', 'beta_save', 'loss_save', 'samp_save', 'betO_save');   % current results
    else
         save(mat_filename, ...
        'IC', 'n', 'p', 'cov', ...  % problem parameters
        'rr', 'rep', 'bb', 'B', ... % current progress indicators
        'rand_save', ... % 'rnd_strm_start', 'rnd_state_last', ... % current rand gen state
        'beta_supp_est_save', 'beta_save', 'loss_save', 'samp_save');   % current results
    end
end

diary off;

end

%% function: support estimation, which does bisection search for the regularization parameter
% Output:
%    beta_ret: beta to return
%    RegPar_ret: regularization parameter to return
%    phi_ret: phi to return
function [phi_ret, beta_ret, RegPar_ret] = support_estimation(L, R, Z, target_nnz, opts, regpar)

[n, p] = size(Z);

loss = intvl_cnsr_cvx_loglikelihood_class_1(n, L, R);
optz = intvl_cnsr_cvx_loglikelihood_maximize_class_1(Z, loss, opts);

assert(isfield(regpar,'Tol') && ~isempty(regpar.Tol), 'regpar.Tol must be set');

% search for target_nnz
beta_nnz_suff = p;  % initialize beta_nnz_suff to p
while (regpar.Hgh - regpar.Low)/geomean([regpar.Hgh regpar.Low]) > regpar.Tol
    % set regularization parameter to the geometric mean
    opts.theta = geomean([regpar.Hgh regpar.Low]);

    % print
    if opts.DisplayLevel >= 0.5
        if ~exist('jj','var')
            jj = 1;
            opts.print.head = strcat(opts.print.head, {'  j    high    curr     low  tgt'});
            existing_len = length(opts.print.text);
            opts.print.text = strcat(opts.print.text, sprintf(' %2d  %3.1e %3.1e %3.1e %3d', jj, regpar.Hgh, opts.theta, regpar.Low, target_nnz));
        else
            jj = jj+1;
            opts.print.text = strcat(opts.print.text(1:existing_len), sprintf(' %2d  %3.1e %3.1e %3.1e %3d', jj, regpar.Hgh, opts.theta, regpar.Low, target_nnz));
        end        
    end

    % call solver for support estimation
    optz = update_regularization(optz, opts);
    optz = update_print(optz, opts);
    [phi, beta] = maximize_by_alt_greedyBCGD(optz);
    beta_nnz = nnz(beta);
    RegPar_ret = opts.theta;

    if beta_nnz > target_nnz
        % increase lower bound
        regpar.Low = opts.theta;
        if beta_nnz < beta_nnz_suff
            % remember this best sufficient beta so far
            beta_ret = beta;
            phi_ret = phi;
            beta_nnz_suff = beta_nnz;
        end
    elseif beta_nnz < target_nnz
        % decrease upper bound
        regpar.Hgh = opts.theta;
    else
        beta_ret = beta;
        phi_ret = phi;
        break;
    end
end

clear loss optz
end

%% function restricted estimation:
function [beta, val, phi] = restricted_estimation(L, R, Z, supp, opts, j_set, initial_beta_S, initial_phi)

    %% support restricted estimation options
    if nargin >= 7 && ~isempty(initial_beta_S)
        opts.Initial_beta = [initial_beta_S;zeros(numel(j_set),1)];
        if nargin >=8 && ~isempty(initial_phi)
            opts.Initial_phi = initial_phi;
        end
    end

    loss = intvl_cnsr_cvx_loglikelihood_class_1(size(Z,1), L, R);
    optz = intvl_cnsr_cvx_loglikelihood_maximize_class_1(Z(:,supp), loss, opts); 
    
    try
        %% call solver
        if nargout == 3
            [val, phi, beta] = maximize_by_knitro(optz, []);
        else
            [val, ~, beta] = maximize_by_knitro(optz, []);
        end

        %% print results
        if isempty(j_set)
            j_set = 0;
        end
        fprintf('%s% 4d % 4d  %3.1e  %3.1e  %3.1e\n', opts.print.text, j_set(1), 0, val, 0, 0);
    catch ME
        clear mex;
        
        val = inf;
        beta = inf(numel(supp),1);
        phi = []; 
        
        %% print results
        if isempty(j_set)
            j_set = 0;
        end
        fprintf('%s% 4d % 4d  %3.1e  %3.1e  %3.1e\n', opts.print.text, j_set(1), 0, val, 0, 0);
        fprintf('\nError message:%s\n', ME.message);
        keyboard;
    end
end



%% read and assign options
function opts = ReadOptions(opts, input_opts)
% 
    opts_size = numel(fieldnames(opts));
    
    if isstruct(input_opts)
        f = fieldnames(input_opts);
        for i = 1:length(f)
            if ~isfield(opts,f{i}); warning('Unrecognized feild: %s\n', f{i}); end
            opts.(f{i}) = input_opts.(f{i});
        end
    else
        if rem(length(input_opts)) == 1
            error('Options must come in as pairs like {''use_parfor'', true, ''estimation_warmstart'', false}')
        end

        for ii = 1:2:length(input_opts)
            if isstring(input_opts(ii+1))
                opts.(input_opts(ii)) = lower(input_opts(ii+1));
            else
                opts.(input_opts(ii)) = input_opts(ii+1);
            end
        end

        if (numel(fieldnames(opts)) > opts_size)
            error('There are unrecognized options.')
        end
    end
end