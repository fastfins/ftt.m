classdef FTTamen < FTT
    %Construct function tensor train using AMEN
            
    properties
        res_x
        res_w
    end
    
    methods
        function obj = FTTamen(func, d, arg, varargin)
            %Construct function tensor train using AMEN
            %        
            obj@FTT(func, d, arg, varargin{:});
            %
            obj.name = 'amen';
            %
            % initialise AMEN
            if size(obj.sample_x, 2) < obj.opt.kick_rank
                disp('Not enough number of user provided samples to enrich ftt')
                obj.sample_x = [];
                obj.res_x{d} = sample_domain(obj.oneds{d}, obj.opt.kick_rank);
                for k = (d-1):-1:1
                    obj.res_x{k} = [obj.res_x{k+1}; sample_domain(obj.oneds{k}, obj.opt.kick_rank)];
                end
            else
                % nested interpolation points for res
                for k = d:-1:1
                    obj.res_x{k} = obj.sample_x(k:d, 1:obj.opt.kick_rank);
                end
                obj.sample_x(:,1:obj.opt.kick_rank) = [];
            end
            % res weights, used on the right
            for k = 1:d
                obj.res_w{k} = randn(obj.opt.init_rank, obj.opt.kick_rank);
            end
            %
            % get dimensions of each TT block
            [d,rs,~] = size(obj);
            % start
            f_evals  = 0;
            als_iter = 0;
            errs     = ones(1,d);
            %
            while true % run ALS
                if obj.direction > 0
                    ind = 1:(d-1);
                else
                    ind = d:-1:2;
                end
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
                        [obj.cores{k}, obj.interp_x{k}, obj.res_w{k}, obj.res_x{k}, obj.cores{k+1}] = build_basis(obj, ...
                            obj.oneds{k}, Jx_left, Jr_left, Rw_left, obj.res_w{k+1}, obj.cores{k+1}, F, Fu, Fr);
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
                        [obj.cores{k}, obj.interp_x{k}, obj.res_w{k}, obj.res_x{k}, obj.cores{k-1}] = build_basis(obj,...
                            obj.oneds{k}, Jx_right, Jr_right, obj.res_w{k-1}, Rw_right, obj.cores{k-1}, F, Fu, Fr);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [core, interp_x, res_w, res_x, core_next] = build_basis(obj, oned, ...
                interp_xold, res_xold, res_w_l, res_w_r, core_next, F, Fu, Fr)
            % find the eigen function of the Schmidt operator
            % cross case, nileft x niright kernels,
            %
            nbleft  = size(F, 1);
            nbright = size(F, 3);
            nrleft  = size(Fr,1);
            nrright = size(Fr,3);
            nnodes  = oned.num_nodes;
            m       = size(F, 4);
            % dimension of the next core
            rn1 = size(core_next, 1);
            nn  = size(core_next, 2);
            rn2 = size(core_next, 3);
            %
            if obj.direction > 0
                % from left >k is integrated
                % T*T' give the schmidt operator evaluated at the quadrature points
                F       = reshape(permute( reshape(F,  nbleft, nnodes, nbright, []), [2,1,3,4]), nnodes*nbleft, []);
                Fu      = reshape(permute( reshape(Fu, nbleft, nnodes, nrright, []), [2,1,3,4]), nnodes*nbleft, []);
                rold    = nbleft;
                rnext   = rn1;
            else
                % <k is integrated
                % T'*T give the schmidt operator evaluated at the quadrature points
                F       = reshape(permute( reshape(F,  nbleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []);
                Fu      = reshape(permute( reshape(Fu, nrleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []);
                rold    = nbright;
                rnext   = rn2;
            end
            %
            [B,A,r] = FTT.local_truncate(obj.opt.loc_err_tol, 1, max(obj.opt.max_rank-obj.opt.kick_rank,1), oned, F);
            %
            % compute residual
            % BSV' = FT is arranged as nodes, rold (fixed), rnew (to be updated)
            % rnew dimension needs to be projected onto res basis from the update
            % direction, rold dimension needs to be projected onto res basis from
            % the fixed dimension
            % for the right projection
            if obj.direction > 0
                % A is r x rnext x m, multiply the rnext dim with res_w_r
                tmp_r   = reshape(permute(reshape(A,r,rn1,m), [1,3,2]), [], rn1)*res_w_r;
                tmp_r   = reshape(permute(reshape(tmp_r,r,m,[]),[1,3,2]),r,[]);
                Fu      = Fu - B*tmp_r;
                % for the left projection
                tmp_l   = res_w_l*reshape(permute(reshape(B,nnodes,rold,r), [2,1,3]), rold, []);
                tmp_l   = reshape(tmp_l, nrleft*nnodes, r);
                % align Fr as rold (nrleft), nodes, rnew (nrright), m
                Fr      = reshape(Fr,nrleft,nnodes,nrright,[]) - reshape(tmp_l*tmp_r,nrleft,nnodes,nrright,[]);
                Fr      = reshape(permute(Fr,[2,1,3,4]), nnodes*nrleft, []);
                rrold   = nrleft;
            else
                % A is r x rnext x m, multiply the rnext dim with res_w_l
                tmp_lt  = res_w_l*reshape(permute(reshape(A,r,rn2,m), [2,1,3]), rn2, []);
                % r x nr_l x m
                tmp_lt  = reshape(permute(reshape(tmp_lt,[],r,m),[2,1,3]),r,[]);
                % Fu is aligned with F
                Fu      = Fu - B*tmp_lt;
                % for the right projection
                tmp_r   = reshape(permute(reshape(B,nnodes,rold,r), [3,1,2]), [], rold)*res_w_r; % r x n x nr_r
                tmp_r   = reshape(tmp_r, r, nnodes*nrright);
                % align Fr as rnew (nrleft), nodes, rold (nrright), m
                Fr      = reshape(Fr,nrleft,nnodes,nrright,m) - permute(reshape(tmp_lt'*tmp_r,nrleft,m,nnodes,nrright), [1,3,4,2]);
                Fr      = reshape(permute(Fr,[2,3,1,4]), nnodes*nrright, []);
                rrold   = nrright;
            end
            % enrich basis
            if isa(oned, 'piecewise')
                % T = kron(speye(rold), oned.mass_L') * [B, Fu];
                T = reshape(oned.mass_L'*reshape([B,Fu],oned.num_nodes,[]),size(B,1),[]);
                [Q,R] = qr(T,0);
                if size(T,1) < size(T,2)
                    Q = Q(:,1:size(T,1));
                end
                % interpolation basis
                % B = full(kron(speye(rold),oned.mass_L') \ Q);
                B = reshape(oned.mass_L'\reshape(Q,oned.num_nodes,[]),size(Q,1),[]);
            else
                %{
                T = [B,Fu];
                [B,R] = qr(T,0);
                if size(T,1) < size(T,2)
                    B = B(:,1:size(T,1));
                end
                %}
                T = reshape(oned.node2basis*reshape([B,Fu],oned.num_nodes,[]),size(B,1),[]);
                [Q,R] = qr(T,0);
                if size(T,1) < size(T,2)
                    Q = Q(:,1:size(T,1));
                end
                B = reshape(oned.basis2node*reshape(Q,oned.num_nodes,[]),size(Q,1),[]);
            end
            rnew        = size(B,2);
            %
            %[indf,core,interp_atx] = point_selection(oned, obj.opt.int_method, reshape(B,nnodes*rold,[]));
            [indf,core,interp_atx] = d_point_selection(obj.opt.int_method, reshape(B,nnodes*rold,[]));
            couple      = reshape(interp_atx*(R(1:rnew,1:r)*A), rnew, rnext, m);
            interp_x    = FTT.local_index(oned, obj.direction, interp_xold, indf);
            %
            Qr      = FTT.local_truncate(obj.opt.loc_err_tol*eps, 1, obj.opt.kick_rank, oned, Fr);
            % rold is the leading dimension in the point selection
            %indr    = point_selection(oned, obj.opt.int_method, reshape(Qr,nnodes*rrold,[]));
            indr    = d_point_selection(obj.opt.int_method, reshape(Qr,nnodes*rrold,[]));
            res_x   = FTT.local_index(oned, obj.direction, res_xold, indr);
            %
            if obj.direction > 0
                % from left >k is integrated
                % get the index for transpose A(:,:,j) for each j
                core = permute(reshape(core,nnodes,rold,rnew), [2,1,3]);
                %the right factor is r x (rn1+enrich) x m, only push the r x rn1 x m
                %block to the next block, first permute the block to m x r x rn1
                core_next = reshape(permute(couple(:,1:rnext,:), [3,1,2]), [], rnext) * reshape(core_next, rnext, []);
                core_next = permute(reshape(core_next, m, rnew, nn, rn2), [2, 3, 4, 1]);
                %
                tmp     = reshape(permute(reshape(res_w_l*reshape(core, rold, nnodes*rnew),rrold,nnodes,rnew),[2,1,3]),[],rnew);
                %if isa(oned, 'spectral')
                %    tmp = kron(speye(rrold),oned.basis2node)*tmp;
                %end
                res_w   = tmp(indr,:);
            else
                % <k is integrated
                % unfold the 3D tensor
                core    = permute(reshape(core,nnodes,rold,rnew), [3,1,2]);
                %the right factor is r x (rn2+enrich) x m, only push the r x rn2 x m
                %block to the next block, first permute the block to rn2 x r x m
                core_next = reshape(core_next, [], rnext) * reshape(permute(couple(:,1:rnext,:), [2,1,3]), rnext, []);
                core_next = reshape(core_next, rn1, nn, rnew, m);
                %
                tmp     = reshape(reshape(core, [], rold)*res_w_r, rnew,[]);
                %if isa(oned, 'spectral')
                %    tmp = tmp*kron(speye(rrold),oned.basis2node)';
                %end
                res_w   = tmp(:,indr);
            end 
        end
    end
end