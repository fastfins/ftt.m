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
                        [obj.cores{k}, obj.interp_x{k}, obj.cores{k+1}] = build_basis(obj, ...
                            obj.oneds{k}, Jx_left, obj.cores{k+1}, cat(3,F,Fe));
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
                        [obj.cores{k}, obj.interp_x{k}, obj.cores{k-1}] = build_basis(obj, ...
                            obj.oneds{k}, Jx_right, obj.cores{k-1}, cat(1,F,Fe));
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
                    disp(rs)
                    break;
                else
                    % flip direction
                    obj.direction = -obj.direction;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [core, interp_x, core_next] = build_basis(obj, oned, interp_xold, core_next, F)
            % find the eigen function of the Schmidt operator
            % cross case, nileft x niright kernels,
            % dimension of the refinement block
            nbleft  = size(F, 1);
            nnodes  = size(F, 2);
            nbright = size(F, 3);
            m       = size(F, 4);
            % dimension of the next core
            rn1 = size(core_next, 1);
            nn  = size(core_next, 2);
            rn2 = size(core_next, 3);
            %
            if obj.direction > 0
                % from left >k is integrated
                % T*T' give the schmidt operator evaluated at the quadrature points
                F       = reshape(permute( reshape(F, nbleft, nnodes, nbright, []), [2,1,3,4]), nnodes*nbleft,  []); %(nright+enrich)*m
                rold    = nbleft;
                %rnext   = rn1;
            else
                % <k is integrated
                % T'*T give the schmidt operator evaluated at the quadrature points
                F       = reshape(permute( reshape(F, nbleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []); %(nbleft+enrich)*m
                rold    = nbright;
                %rnext   = rn2;
            end
            %
            [B,A,r]     = FTT.local_truncate(obj.opt.loc_err_tol, 1, obj.opt.max_rank, oned, F);
            %[ind,core,interp_atx] = point_selection(oned, obj.opt.int_method, reshape(B,nnodes*rold,[]));
            [ind,core,interp_atx] = d_point_selection(obj.opt.int_method, reshape(B,nnodes*rold,[]));
            couple      = reshape(interp_atx*A, r, [], m);
            interp_x    = FTT.local_index(oned, obj.direction, interp_xold, ind);
            %
            if obj.direction > 0
                % from left >k is integrated
                % get the index for transpose A(:,:,j) for each j
                core = permute(reshape(core,nnodes,rold,r), [2,1,3]);
                %the right factor is r x (rn1+enrich) x m, only push the r x rn1 x m
                %block to the next block, first permute the block to m x r x rn1
                core_next = reshape(permute(couple(:,1:rn1,:), [3,1,2]), [], rn1) * reshape(core_next, rn1, []);
                core_next = permute(reshape(core_next, m, r, nn, rn2), [2, 3, 4, 1]);
            else
                % <k is integrated
                % unfold the 3D tensor
                core    = permute(reshape(core,nnodes,rold,r), [3,1,2]);
                %the right factor is r x (rn2+enrich) x m, only push the r x rn2 x m
                %block to the next block, first permute the block to rn2 x r x m
                core_next = reshape(core_next, [], rn2) * reshape(permute(couple(:,1:rn2,:), [2,1,3]), rn2, []);
                core_next = reshape(core_next, rn1, nn, r, m);
            end
        end
    end
end