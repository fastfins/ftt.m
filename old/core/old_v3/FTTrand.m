classdef FTTrand < FTT
    %Construct function tensor train using random enrichment
    
    methods
        function obj = FTTrand(func, d, arg, varargin)
            %Construct function tensor train using random enrichment
            %
            obj@FTT(func, d, arg, varargin{:});
            %
            obj.name = 'rand';
            %
            % get dimensions of each TT block
            [d,rs,~] = size(obj);
            % start
            f_evals  = 0;
            als_iter = 0;
            errs     = ones(1,d);
            while true % run ALS
                % at the head, update the random enrichment set
                if size(obj.sample_x, 2) < obj.opt.kick_rank
                    obj.sample_x = [];
                    enrich  = zeros(d, obj.opt.kick_rank);
                    for k = 1:d
                        enrich(k,:) = sample_domain(obj.oneds{k}, obj.opt.kick_rank);
                    end
                else
                    enrich = obj.sample_x(:,1:obj.opt.kick_rank);
                    obj.sample_x(:,1:obj.opt.kick_rank) = [];
                end
                %
                if obj.direction > 0
                    ind = 1:(d-1);
                else
                    ind = d:-1:2;
                end
                %
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
                            obj.direction, obj.opt.int_method, obj.opt.loc_err_tol, obj.opt.max_rank);
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
                            obj.direction, obj.opt.int_method, obj.opt.loc_err_tol, obj.opt.max_rank);
                        rs(k-1) = size(obj.cores{k}, 1);
                    end
                end
                als_iter = als_iter + 1;
                % evaluate the ftt and give error estimates
                debug_errs = [];
                if ~isempty(obj.debug_x)
                    exact   = func(obj.debug_x);
                    approx  = eval(obj, obj.debug_x);
                    debug_errs(1) = max(abs(exact(:) - approx(:)))  / max(abs(exact(:)));
                    debug_errs(2) = mean(abs(exact(:) - approx(:))) / max(abs(exact(:)));
                end
                FTT.print_iteration(als_iter, errs, rs, f_evals, debug_errs);
                %
                if als_iter == obj.opt.max_als || max(errs) < obj.opt.err_tol
                    disp('ALS completed')
                    break;
                else
                    % flip direction
                    obj.direction = -obj.direction;
                end
            end
            obj = round(obj);
        end
    end
end