function obj = cross(obj, func, d, sample_x, debug_x)
% cross iterations

% initilise TT cores
if isempty(obj.cores)
    % copy properties
    obj.direction = 1;
    % nested interpolation points
    obj.cores = cell(d,1);
    %
    % use the user provided sample set
    if size(sample_x, 2) < obj.opt.init_rank
        disp('Not enough number of samples to initialise ftt')
        sample_x = [];
        obj.interp_x{d} = sample_domain(obj.oneds{d}, obj.opt.init_rank);
        for k = (d-1):-1:1
            obj.interp_x{k} = [obj.interp_x{k+1}; sample_domain(obj.oneds{k}, obj.opt.init_rank)];
        end
    else
        for k = d:-1:1
            obj.interp_x{k} = sample_x(k:d, 1:obj.opt.init_rank);
        end
        sample_x(:,1:obj.opt.init_rank) = [];
    end
    %
    % determine m
    y = func(obj.interp_x{k}(:,1));
    m = size(y,1);
    % interpolation weights
    obj.cores{1} = zeros(1, obj.oneds{1}.num_nodes, obj.opt.init_rank, m);
    for k = 2:d-1
        obj.cores{k} = zeros(obj.opt.init_rank, obj.oneds{k}.num_nodes, obj.opt.init_rank);
    end
    obj.cores{d} = zeros(obj.opt.init_rank, obj.oneds{d}.num_nodes, 1);
    % initialise the residual blocks for AMEN
    if strcmp(obj.opt.tt_method, 'amen')
        if size(sample_x, 2) < obj.opt.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            obj.res_x{d} = sample_domain(obj.oneds{d}, obj.opt.kick_rank);
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
else
    obj.direction = -obj.direction;
end
% reinitialise the residual blocks for AMEN
if strcmp(obj.opt.tt_method, 'amen') && isempty(obj.res_w)
    if obj.direction < 0 % direction has already been flipped
        for k = 1:(d-1)
            obj.res_w{k} = randn(obj.opt.kick_rank, size(obj.cores{k},3));
        end
        obj.res_w{d} = randn(size(obj.cores{d},1), obj.opt.kick_rank);
    else
        obj.res_w{1} = randn(obj.opt.kick_rank, size(obj.cores{1},3));
        for k = 2:d
            obj.res_w{k} = randn(size(obj.cores{k},1), obj.opt.kick_rank);
        end
    end
end

% get dimensions of each TT block
% start
f_evals  = 0;
als_iter = 0;
errs = zeros(1,d);
rs = ones(1,d);
while true % run ALS
    if obj.direction > 0
        ind = 1:(d-1);
    else
        ind = d:-1:2;
    end
    %
    switch obj.opt.tt_method
        case {'random'}
            % at the head, update the random enrichment set
            if size(sample_x, 2) < obj.opt.kick_rank
                sample_x = [];
                enrich  = zeros(d, obj.opt.kick_rank);
                for k = 1:d
                    enrich(k,:) = sample_domain(obj.oneds{k}, obj.opt.kick_rank);
                end
            else
                enrich = sample_x(:,1:obj.opt.kick_rank);
                sample_x(:,1:obj.opt.kick_rank) = [];
            end
            % start
            for k = ind
                if obj.direction > 0
                    if k == 1
                        Jx_left = [];
                    else
                        Jx_left = obj.interp_x{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, Jx_left, obj.interp_x{k+1}, func);
                    [Fe,ne] = FTT.local_block(obj.opt.sqrt_flag, obj.oneds{k}, Jx_left, enrich(k+1:d,:),   func);
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
                    if k == d
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
                    if k == d
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
                    if k < d
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
        debug_errs(1) = max(abs(exact(:) - approx(:)))  / max(abs(exact(:)));
        debug_errs(2) = mean(abs(exact(:) - approx(:))) / max(abs(exact(:)));
    end
    % Print the information of each TT iteration
    fprintf('als=%2d, max_local_error=%3.3e, mean_local_error=%3.3e, max_rank=%d, cum#fevals=%3.3e\n', ...
        als_iter, max(errs), mean(errs), max(rs), f_evals);
    if ~isempty(debug_errs)
        fprintf('als=%2d, max_debug_error=%3.3e, mean_debug_error=%3.3e\n', ...
            als_iter, debug_errs(1), debug_errs(2));
    end
    fprintf('\n');
    %
    if als_iter == obj.opt.max_als || max(errs) < obj.opt.als_tol
        disp('ALS completed')
        break;
    else
        % flip direction
        obj.direction = -obj.direction;
    end
end
disp(rs)
end