function obj = cross(obj, func, sample_x, debug_x)
% Cross iterations for building the TT, only call this function if you know
% what you are doing

% initilise TT cores
if isempty(obj.cores)
    % copy properties
    obj.direction = 1;
    % nested interpolation points
    obj.cores = cell(obj.d,1);
    %
    % use the user provided sample set
    if size(sample_x, 2) < obj.opt.init_rank
        disp('Not enough number of samples to initialise ftt')
        sample_x = [];
        obj.interp_x{obj.d} = sample_domain(obj.oneds{obj.d}, obj.opt.init_rank);
        for k = (obj.d-1):-1:1
            obj.interp_x{k} = [obj.interp_x{k+1}; sample_domain(obj.oneds{k}, obj.opt.init_rank)];
        end
    else
        for k = obj.d:-1:1
            obj.interp_x{k} = sample_x(k:obj.d, 1:obj.opt.init_rank);
        end
        sample_x(:,1:obj.opt.init_rank) = [];
    end
    %
    % determine m
    y = func(obj.interp_x{k}(:,1));
    m = size(y,1);
    % interpolation weights
    obj.cores{1} = zeros(1, obj.oneds{1}.num_nodes, obj.opt.init_rank, m);
    for k = 2:obj.d-1
        obj.cores{k} = zeros(obj.opt.init_rank, obj.oneds{k}.num_nodes, obj.opt.init_rank);
    end
    obj.cores{obj.d} = zeros(obj.opt.init_rank, obj.oneds{obj.d}.num_nodes, 1);
    %{
    % initialise the residual blocks for AMEN
    if strcmp(obj.opt.tt_method, 'amen') && isempty(obj.res_x)
        if size(sample_x, 2) < obj.opt.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            obj.res_x{obj.d} = sample_domain(obj.oneds{obj.d}, obj.opt.kick_rank);
            for k = (d-1):-1:1
                obj.res_x{k} = [obj.res_x{k+1}; sample_domain(obj.oneds{k}, obj.opt.kick_rank)];
            end
        else
            % nested interpolation points for res
            for k = d:-1:1
                obj.res_x{k} = sample_x(k:d, 1:obj.opt.kick_rank);
            end
        end
        % res weights, used on the right
        for k = 1:d
            obj.res_w{k} = randn(obj.opt.init_rank, obj.opt.kick_rank);
        end
    end
    %}
else
    obj.direction = -obj.direction;
    %{
    if obj.direction > 0
        if size(sample_x, 2) < obj.opt.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            obj.res_x{obj.d} = sample_domain(obj.oneds{obj.d}, obj.opt.kick_rank);
            for k = (d-1):-1:1
                obj.res_x{k} = [obj.res_x{k+1}; sample_domain(obj.oneds{k}, obj.opt.kick_rank)];
            end
        else
            % nested interpolation points for res
            for k = d:-1:1
                obj.res_x{k} = sample_x(k:d, 1:obj.opt.kick_rank);
            end
        end
    else
    end
    % reinitialise the residual blocks for AMEN
    if strcmp(obj.opt.tt_method, 'amen') && isempty(obj.res_w)
        if obj.direction < 0 % direction has already been flipped
            for k = 1:(d-1)
                obj.res_w{k} = randn(obj.opt.kick_rank, size(obj.cores{k},3));
            end
            obj.res_w{obj.d} = randn(size(obj.cores{obj.d},1), obj.opt.kick_rank);
        else
            obj.res_w{1} = randn(obj.opt.kick_rank, size(obj.cores{1},3));
            for k = 2:d
                obj.res_w{k} = randn(size(obj.cores{k},1), obj.opt.kick_rank);
            end
        end
    end
    %}
end
% initialise the residual coordinates for AMEN
if strcmp(obj.opt.tt_method, 'amen') && isempty(obj.res_x)
    if obj.direction > 0 % direction has already been flipped
        if size(sample_x, 2) < obj.opt.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            obj.res_x{obj.d} = sample_domain(obj.oneds{obj.d}, obj.opt.kick_rank);
            for k = (obj.d-1):-1:1
                obj.res_x{k} = [obj.res_x{k+1}; sample_domain(obj.oneds{k}, obj.opt.kick_rank)];
            end
        else
            % nested interpolation points for res
            for k = obj.d:-1:1
                obj.res_x{k} = sample_x(k:obj.d, 1:obj.opt.kick_rank);
            end
        end
    else
        if size(sample_x, 2) < obj.opt.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            obj.res_x{1} = sample_domain(obj.opt.kick_rank, obj.oneds{1});
            for k = 2:obj.d
                obj.res_x{k} = [sample_domain(obj.opt.kick_rank, obj.oneds{k}); obj.res_x{k-1}];
            end
        else
            % nested interpolation points for res
            for k = obj.d:-1:1
                obj.res_x{k} = sample_x(1:obj.opt.kick_rank, k:obj.d);
            end
        end
    end
end
% reinitialise the residual blocks for AMEN
if strcmp(obj.opt.tt_method, 'amen') && isempty(obj.res_w)
    if obj.direction < 0 % direction has already been flipped
        for k = 1:(obj.d-1)
            obj.res_w{k} = randn(obj.opt.kick_rank, size(obj.cores{k},3));
        end
        obj.res_w{obj.d} = randn(size(obj.cores{obj.d},1), obj.opt.kick_rank);
    else
        obj.res_w{1} = randn(obj.opt.kick_rank, size(obj.cores{1},3));
        for k = 2:obj.d
            obj.res_w{k} = randn(size(obj.cores{k},1), obj.opt.kick_rank);
        end
    end
end
% 
if isempty(debug_x)
    fprintf('\n\n  >> ALS  max_local_E  mean_local_E  max_r  cum#fevals\n');
else
    fprintf('\n\n  >> ALS  max_local_E  mean_local_E  max_r  max_debug_E  mean_debug_E  cum#fevals\n');
end

    
% start
l2_err = 0;
f_evals = 0;
als_iter = 0;
rs = ones(1,obj.d);
while true % run ALS
    if obj.direction > 0
        ind = 1:(obj.d-1);
    else
        ind = obj.d:-1:2;
    end
    %
    switch obj.opt.tt_method
        case {'random'}
            % at the head, update the random enrichment set
            if size(sample_x, 2) < obj.opt.kick_rank
                sample_x = [];
                enrich  = zeros(obj.d, obj.opt.kick_rank);
                for k = 1:obj.d
                    enrich(k,:) = sample_domain(obj.oneds{k}, obj.opt.kick_rank);
                end
            else
                enrich = sample_x(:,1:obj.opt.kick_rank);
                sample_x(:,1:obj.opt.kick_rank) = [];
            end
            % start
            errs = zeros(1,obj.d);
            for k = ind
                if obj.direction > 0
                    if k == 1
                        Jx_left = [];
                    else
                        Jx_left = obj.interp_x{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, Jx_left, obj.interp_x{k+1}, func);
                    [Fe,ne] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, Jx_left, enrich(k+1:obj.d,:),   func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    f_evals = f_evals + nf + ne;
                    % Schmidt operator and interpolation, k - 1 is the previous index
                    % push couple matrix to the right
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k+1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_left, obj.cores{k+1}, cat(3,F,Fe), ...
                        obj.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank);
                    rs(k) = size(obj.cores{k}, 3);
                else
                    if k == obj.d
                        Jx_right = [];
                    else
                        Jx_right = obj.interp_x{k+1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, obj.interp_x{k-1}, Jx_right, func);
                    % different from left iteration, k + 1 is the previous index
                    [Fe,ne] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, enrich(1:k-1,:),   Jx_right, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    f_evals = f_evals + nf + ne;
                    % Schmidt operator and interpolation, k + 1 is the previous index
                    % push couple matrix to the right
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k-1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_right, obj.cores{k-1}, cat(1,F,Fe), ...
                        obj.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank);
                    rs(k-1) = size(obj.cores{k}, 1);
                end
            end
        case {'amen'}
            % start
            errs = zeros(1,obj.d);
            for k = ind
                if obj.direction > 0
                    if k == 1
                        Jx_left = [];
                        Jr_left = [];
                        Rw_left = 1;
                    else
                        Jx_left = obj.interp_x{k-1};
                        Jr_left = obj.res_x{k-1};
                        Rw_left = obj.res_w{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, Jx_left, obj.interp_x{k+1}, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    % evaluate res function at x_k nodes
                    [Fr,nr] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, Jr_left, obj.res_x{k+1}, func);
                    f_evals = f_evals + nf + nr;
                    % evaluate update function at x_k nodes
                    if k > 1
                        [Fu,nu] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, Jx_left, obj.res_x{k+1}, func);
                        f_evals = f_evals + nu;
                    else
                        Fu = Fr;
                    end
                    [obj.cores{k}, obj.interp_x{k}, obj.res_w{k}, obj.res_x{k}, obj.cores{k+1}] = FTT.build_basis_amen(...
                        obj.oneds{k}, Jx_left, Jr_left, Rw_left, obj.res_w{k+1}, obj.cores{k+1}, F, Fu, Fr, ...
                        obj.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank, obj.opt.kick_rank);
                    rs(k) = size(obj.cores{k}, 3);
                else
                    if k == obj.d
                        Jx_right = [];
                        Jr_right = [];
                        Rw_right = 1;
                    else
                        Jx_right = obj.interp_x{k+1};
                        Jr_right = obj.res_x{k+1};
                        Rw_right = obj.res_w{k+1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, obj.interp_x{k-1}, Jx_right, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    % evaluate res function at x_k nodes
                    [Fr,nr] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, obj.res_x{k-1}, Jr_right, func);
                    f_evals = f_evals + nf + nr;
                    % different from left iteration, k + 1 is the previous index
                    % evaluate update function at x_k nodes
                    if k < obj.d
                        [Fu,nu] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, obj.res_x{k-1}, Jx_right, func);
                        f_evals = f_evals + nu;
                    else
                        Fu = Fr;
                    end
                    [obj.cores{k}, obj.interp_x{k}, obj.res_w{k}, obj.res_x{k}, obj.cores{k-1}] = FTT.build_basis_amen(...
                        obj.oneds{k}, Jx_right, Jr_right, obj.res_w{k-1}, Rw_right, obj.cores{k-1}, F, Fu, Fr, ...
                        obj.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank, obj.opt.kick_rank);
                    rs(k-1) = size(obj.cores{k}, 1);
                end
            end
    end
    als_iter = als_iter + 1;
    % evaluate the ftt and give error estimates
    debug_errs = [];
    if ~isempty(debug_x)
        exact   = func(debug_x);
        approx  = eval(obj, debug_x);
        l2_err = mean((exact(:) - approx(:)).^2).^0.5;
        debug_errs(1) = max(abs(exact(:) - approx(:)))  / max(abs(exact(:)));
        debug_errs(2) = l2_err / mean(exact(:).^2).^0.5;
    end
    % Print the information of each TT iteration
    if isempty(debug_errs)
        fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10i\n', ...
            [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), f_evals]);
    else
        fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10.5E   %10.5E  %10i\n', ...
            [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), debug_errs(1), debug_errs(2), f_evals]);
    end
    %
    if als_iter == obj.opt.max_als || max(errs(ind)) < obj.opt.als_tol
        disp('  >> ALS completed, TT ranks')
        disp(['  >> ' num2str(rs)])
        break;
    else
        % flip direction
        obj.direction = -obj.direction;
    end
end
obj.l2_err = l2_err;
obj.n_evals = f_evals;
end