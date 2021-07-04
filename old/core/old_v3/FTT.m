classdef FTT
    % FTT abstract class, need to use either FTTamen or FTTrand
    %
    % FTT Properties:
    %
    %   opt         - FTToption
    %
    %   direction   - ALS direction, >0: built from left to right
    %                                <0: built from right to left
    %
    %   oneds{k}    - data structure containing information for building
    %                 one dimensional polynomial representation
    %
    %   cores{k}    - nodal values or coefficent tensors of oned functions
    %                 the dimension of the current core is organised as
    %                 previous rank x oned{k}.num_nodes x current rank
    %
    %   interp_x{x} - interpolation coordinates
    %
    %   isrounded   - indicate if the TT is rounded after construction
    %
    %
    % FTToption Methods:
    %
    %   eval        - evaluate fTT. The output is horizontally aligned.
    %
    %   eval_block  - evaluate fTT for either the first or last k variables.
    %                 The output is horizontally aligned.
    %
    %   round       - round the TT cores
    %
    %   int         - integrate the entire TT
    %
    %   int_block   - integrate a block of TT cores
    %
    %   size        - size of the TT
    %
    %   get_name    - get the name of the TT construction method
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example (vector function outputs, m = 2):
    %
    % poly = setup_oned(5, 'type', 'Lagrange', 'lag_elems', 10, 'ghost_size', 1E-10);
    % func = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
    % d    = 10;
    %
    % debug_size = 1E4;
    % debug_x = zeros(d, debug_size);
    % for k = 1:d
    %     debug_x(k,:) = sample_oned_domain(poly, debug_size);
    % end
    %
    % % use alternating energy enrichment (AMEN)
    % opt1 = ftt_options('method', 'AMEN', 'oned_ref', poly, ...
    %     'err_tol', 1E-8, 'loc_err_tol', 1E-10, 'max_rank', 50);
    % ftt1 = build_ftt(func, d, [], opt1, 'debug_x', debug_x);
    %
    % % use random enrichment
    % opt2 = ftt_options('method', 'Random', 'oned_ref', poly, ...
    %     'err_tol', 1E-8, 'loc_err_tol', 1E-10, 'max_rank', 50);
    % ftt2 = build_ftt(func, d, [], opt2, 'debug_x', debug_x);
    %
    % % evaluate the function and its factorisations
    % exact   = func(debug_x);
    % approx1 = eval_ftt(ftt1, debug_x);
    % approx2 = eval_ftt(ftt2, debug_x);
    % % plot the error
    % figure
    % plot(exact(:) - approx1(:), 'x')
    % hold on
    % plot(exact(:) - approx2(:), '.')
    %
    %
    
    properties
        opt FTToption
        cores
        oneds
        interp_x
        direction
        name
        sample_x
        debug_x
        res_x
        res_w
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static)
        
        function [f, f_evals] = local_block(sqrt_flag, oned, xleft, xright, func)
            % Evaluate the user function at cross indices and qudrature or
            % intepolation points weighed by the interplation matrix
            %
            if isempty(xleft)
                % left boudnary
                nileft  = 1;
                niright = size(xright, 2);
                % f       = zeros(niright, oned_def.nquad);
                % define parameters, each xright binding with a block of xquad
                param   = [ repmat(oned.nodes(:)', 1, niright); ...
                    reshape(repmat(xright, oned.num_nodes, 1), size(xright,1), niright*oned.num_nodes) ];
                % return of func is a vector, reshape to a matrix
                % 1st index: xquad, 2nd: xright
                f       = reshape( func(param)', 1, oned.num_nodes, niright, []);
            elseif isempty(xright)
                % right boundary
                nileft  = size(xleft, 2);
                niright = 1;
                % f       = zeros(nileft, oned_def.nquad);
                % define parameters, each xquad binding with a block of xleft
                param   = [ repmat(xleft, 1, oned.num_nodes); ...
                    reshape( repmat(oned.nodes(:)', nileft, 1), 1, nileft*oned.num_nodes ) ];
                % return of func is a vector, reshape to a matrix
                % 1st index: xleft, 2nd: xquad
                f       = reshape( func(param)', nileft, oned.num_nodes, 1, []);
            else
                %
                nileft  = size(xleft, 2);
                niright = size(xright, 2);
                % f       = zeros(nileft, niright, oned_def.nquad);
                
                tmp_lq  = [ repmat(xleft, 1, oned.num_nodes); ...
                    reshape( repmat(oned.nodes(:)', nileft, 1), 1, nileft*oned.num_nodes ) ];
                
                param   = [ repmat(tmp_lq, 1, niright); ...
                    reshape(repmat(xright, oned.num_nodes*nileft, 1), size(xright,1), nileft*oned.num_nodes*niright) ];
                
                % return of func is a vector, reshape to tensor
                f       = reshape( func(param)', nileft, oned.num_nodes, niright, []);
            end
            if sqrt_flag
                if any(f(:)<0)
                    error('Error: negative function evaluations')
                end
                f   = f.^(0.5);
            end
            if isa(oned, 'spectral')
                f = oned.node2basis*reshape(permute(f, [2,1,3,4]), oned.num_nodes, []);
                f = permute(reshape(f, oned.num_nodes, nileft, niright, []), [2 1 3 4]);
            end
            f_evals = size(param, 2);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [B,A,r] = local_truncate(loc_err_tol, min_rank, max_rank, oned, F)
            % Truncate the svd for each TT block
            %
            if isa(oned, 'piecewise')
                %[U,S,V] = svd( kron(speye(rold), obj.oneds{k}.mass_L') * F, 0 );
                [U,S,V] = svd( reshape(oned.mass_L'*reshape(F,oned.num_nodes,[]),size(F,1),[]), 0);
                s = diag(S);
                % truncation index r
                ind = s/s(1) > loc_err_tol;
                r = max(min_rank, min(sum(ind), max_rank));
                % interpolation basis
                %B   = full(kron(speye(rold), obj.oneds{k}.mass_L') \ U(:,1:r));
                B = reshape(oned.mass_L'\reshape(U(:,1:r),oned.num_nodes,[]),size(F,1),[]);
                A = s(1:r).*V(:,1:r)';
            else
                [U,S,V] = svd(F, 0);
                s = diag(S);
                % truncation index r
                ind = s/s(1) > loc_err_tol;
                r = max(min_rank, min(sum(ind), max_rank));
                B = U(:,1:r);
                A = s(1:r).*V(:,1:r)';
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function rerr = local_error(core, f)
            % Compute local error for TT block
            %
            diff = core(:) - f(:);
            rerr = max(abs(diff)) / max(abs(f(:)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function interp_x = local_index(oned, direction, interp_xold, ind)
            % Build the local nested index set
            %
            if isempty(interp_xold)
                % interpolation points
                interp_x = reshape( oned.nodes(ind), 1, [] );
            else
                % the basis coefficients is now an unfolded tensor, npoly x rold x rnew
                % the candidate index set ( x_old, nodes ), nodes is the leading dimension
                nnodes  = oned.num_nodes;
                rold    = size(interp_xold, 2);
                ipair   = [ reshape( repmat(1:rold, nnodes, 1), 1, nnodes*rold); repmat(1:nnodes, 1, rold) ];
                iselect = ipair(:, ind);
                if direction > 0
                    interp_x = [ interp_xold(:, iselect(1,:)); reshape( oned.nodes(iselect(2,:)), 1, length(ind)) ];
                else
                    interp_x = [ reshape( oned.nodes(iselect(2,:)), 1, length(ind)); interp_xold(:, iselect(1,:)) ];
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function print_iteration(iter, local_err, rs, f_evals, debug_errs)
            % Print the information of each TT iteration
            %
            fprintf('als=%2d, max_local_error=%3.3e, mean_local_error=%3.3e, max_rank=%d, cum#fevals=%3.3e\n', ...
                iter, max(local_err), mean(local_err), max(rs), f_evals);
            if ~isempty(debug_errs)
                fprintf('als=%2d, max_debug_error=%3.3e, mean_debug_error=%3.3e\n', ...
                    iter, debug_errs(1), debug_errs(2));
            end
            fprintf('\n');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function obj = FTT(func, d, arg, varargin)
            % Construct tensor train for a function mapping from R^d to R^m.
            %
            % Inputs:
            %
            %   func - a function (R^d to R^m) that take inputs as a dxn
            %          matrix and returns mxn vector
            %
            %   d    - dimension of the input variable
            %
            %   arg  - either an existing ftt used as the initial guess, or
            %          a set of one dimensional bases for discretising the
            %          function. If only one set of basis is supplied, each
            %          coordinate is discretised using the same basis
            %
            %   opt  - FTT options
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % parsing inputs
            defaultOption    = FTToption();
            defaultSampleSet = [];
            defaultDebugSet  = [];
            
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            %
            addRequired(p,'func',@(x) isa(x, 'function_handle'));
            addRequired(p,'d',validScalarPosNum);
            addRequired(p,'arg');
            addOptional(p,'option',defaultOption);
            %
            addParameter(p,'sample_x', defaultSampleSet);
            addParameter(p,'debug_x',  defaultDebugSet);
            %
            p.KeepUnmatched = false;
            parse(p,func,d,arg,varargin{:});
            %
            obj.opt = p.Results.option;
            sample_x    = p.Results.sample_x;
            debug_x     = p.Results.debug_x;
            %
            has_ftt = true;
            if isa(arg, 'cell')
                obj.oneds = arg;
                has_ftt = false;
            else
                if isa(arg, 'FTT')
                    obj = arg;
                    obj.direction = -obj.direction;
                elseif isa(arg, 'oned')
                    for k = 1:d
                        obj.oneds{k} = arg;
                    end
                    has_ftt = false;
                else
                    error('wrong type of argument')
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % initilise
            if ~has_ftt
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
                %
                % initialise the residul blocks for AMEN
                if strcomp (obj.opt.method, 'amen')
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
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cross iterations
            % get dimensions of each TT block
            [~,rs,~] = size(obj);
            % start
            f_evals  = 0;
            als_iter = 0;
            errs     = ones(1,d);
            while true % run ALS
                if obj.direction > 0
                    ind = 1:(d-1);
                else
                    ind = d:-1:2;
                end
                %
                switch obj.opt.method
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
                                    obj.direction, obj.opt.int_method, obj.opt.loc_err_tol, obj.opt.max_rank);
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
                                    obj.direction, obj.opt.int_method, obj.opt.loc_err_tol, obj.opt.max_rank);
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % round the tensor train for amen
            if strcomp (obj.opt.method, 'amen')
                obj = round(obj);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = round(obj, thres)
            % round the TT cores
            %
            if nargin == 1
                thres = obj.opt.loc_err_tol;
            end
            [d,rs,~] = size(obj);
            %
            obj.direction = -obj.direction;
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
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k+1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_left, obj.cores{k+1}, obj.cores{k}, ...
                        obj.direction, obj.opt.int_method, thres, obj.opt.max_rank);
                    rs(k) = size(obj.cores{k}, 3);
                else
                    if k == d
                        Jx_right = [];
                    else
                        Jx_right = obj.interp_x{k+1};
                    end
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k-1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_right, obj.cores{k-1}, obj.cores{k}, ...
                        obj.direction, obj.opt.int_method, thres, obj.opt.max_rank);
                    rs(k-1) = size(obj.cores{k}, 1);
                end
            end
            disp(rs)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = intgrate(obj)
            % Integrate the entire TT
            %
            if obj.direction > 0
                % last dimenion
                figeqk  = 1;
                for k = d:-1:1
                    nx  = obj.oneds{k}.num_nodes;
                    rkm = size(obj.cores{k}, 1);
                    rk  = size(obj.cores{k}, 3);
                    % push the previous dimenion into the current ftt by modifying the
                    % coefficients, then ys{k} is used to replace the ftt core at the
                    % k-th coordinate for obtaining the integrated function over the
                    % last d-k dimensions, ys{k} is a rk-1 by nx matrix
                    ys = reshape(reshape(obj.cores{k}, rkm*nx, rk)*figeqk, rkm, nx);
                    figeqk  = oned_integral(obj.oneds{k}, ys')';
                end
                z = figeqk;
            else
                % first dimenion
                fileqk = 1;
                for k = 1:d
                    nx  = obj.oneds{k}.num_nodes;
                    rkm = size(obj.cores{k}, 1);
                    rk  = size(obj.cores{k}, 3);
                    % push the previous dimenion into the current ftt
                    % ys{k} is a nx by rk matrix
                    ys = reshape(fileqk*reshape(obj.cores{k}, rkm, nx*rk), nx, rk);
                    fileqk  = oned_integral(obj.oneds{k}, ys);
                end
                z = fileqk;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ftt = integrate_block(obj, ind)
            % Integrate a block of TT cores
            %
            %   ind - indices of coordinates will be integrated
            %
            %   ftt - output TT after integration
            %
            ftt = obj;
            d = length(ftt.cores);
            %
            if length(ind) ~= length(unique(ind))
                disp('integration indices are not unique')
                ind = unique(ind);
            end
            %
            if max(ind) > d || min(ind) < 1
                error('integration indices out of bounds')
            end
            %
            if length(ind) == d
                error('integration over all indices, should use int for integration')
            end
            ind = reshape(ind, 1, []);
            if ind(end) == d
                jr = d;
                for i = (length(ind)-1):-1:1
                    j = ind(i);
                    if j+1 == jr
                        % if the index is continous, move to the left
                        jr = j;
                    else
                        % otherwise stop searching, jr is the last stopping index
                        % jr:d are continous block
                        break;
                    end
                end
                ind2 = jr:d;
                ind1 = setdiff(ind, ind2);
                ind2 = fliplr(ind2);
            else
                ind1 = ind;
                ind2 = [];
            end
            %the left blocks
            if ~isempty(ind1)
                for i = 1:length(ind1)
                    k   = ind1(i);
                    nx  = ftt.oneds{k}.num_nodes;
                    
                    rkm = size(ftt.cores{k}, 1);
                    rk  = size(ftt.cores{k}, 3);
                    tmp = oned_integral(ftt.oneds{k}, reshape(permute(ftt.cores{k}, [2,1,3]), nx,rkm*rk));
                    tmp = reshape(tmp, rkm, rk);
                    ftt.cores{k+1} = reshape(tmp*reshape(ftt.cores{k+1}, rk, []), rkm, nx, []);
                end
            end
            % the right blocks
            if ~isempty(ind2)
                for i = 1:length(ind2)
                    k = ind2(i);
                    nx  = ftt.oneds{k}.num_nodes;
                    rkm = size(ftt.cores{k}, 1);
                    rk  = size(ftt.cores{k}, 3); % rk should always be 1
                    tmp = oned_integral(ftt.oneds{k}, reshape(permute(ftt.cores{k}, [2,1,3]), nx,rkm*rk));
                    tmp = reshape(tmp, rkm, rk);
                    ftt.cores{k-1} = reshape(reshape(ftt.cores{k-1}, [], rkm)*tmp, [], nx, rk);
                end
            end
            % delete marginalised blocks
            ftt.oneds(ind) = [];
            ftt.cores(ind) = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fx = eval(obj, x, varargin)
            % Evaluate the ftt function. The output is horizontally aligned.
            %
            %   x  - input variables, d x n.
            %
            %   pt - (optional) permuation index, where fTT is built for f(y),
            %        y is the permutation of x s.t. x = y(p) and y = x(pt)
            %        for a given p, pt can be obtained by pt(p) = 1:length(p);
            %
            %   fx - function values at x, mxn
            %
            
            % permute the input
            if length(varargin) > 1
                error('too many optional arguments')
            elseif length(varargin) == 1
                pt = varargin{1};
                x = x(pt,:);
            end
            fx  = eval_block(obj, x);
            if obj.direction > 0
                fx  = fx';
            end
            if obj.opt.sqrt_flag
                fx  = fx.^2;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fx = eval_block(obj, x)
            % Evaluate the fTT for either the first or last k variables.
            % The output is horizontally aligned.
            %
            %   x  - input variables, dxn.
            %
            %   fx - function values at x, mxn
            %
            d   = length(obj.cores);
            k   = size(x, 1);
            nx  = size(x, 2);
            %
            if obj.direction > 0
                % start from the 1st dimension
                fx  = ones(nx,1);
                % all the intermediate dimensions, except the last dimension
                for j = 1:min(k,d)
                    nj  = size(obj.cores{j}, 2);
                    rjm = size(obj.cores{j}, 1);
                    %
                    if j < d || (size(obj.cores{j}, 3) > 1 && size(obj.cores{j}, 4) == 1)
                        tmp = reshape(permute(obj.cores{j}, [2,1,3]), nj, []);
                    else
                        % collapse the third dimension = 1
                        tmp = reshape(permute(reshape(obj.cores{j}, rjm, nj, []), [2,1,3]), nj, []);
                    end
                    T   = reshape(permute(reshape( eval(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
                    % how to speed up this part?
                    jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
                    ii  = repmat((1:nx)', 1, rjm);
                    B   = sparse(ii(:), jj(:), fx(:), nx, rjm*nx);
                    %
                    fx  = B*T;
                end
            else
                % start from the last dimension
                xind = k:-1:1;
                tind = d:-1:1;
                fx   = ones(1,nx);
                % all the intermediate dimensions, except the first dimension
                % need to walk through d-k+1 dimensions
                for i = 1:min(k,d)
                    j = tind(i);
                    nj  = size(obj.cores{j}, 2);
                    rj  = size(obj.cores{j}, 3);
                    %
                    if j > 1 || (size(obj.cores{j}, 1) > 1 && size(obj.cores{j}, 4) == 1)
                        tmp = reshape(permute(obj.cores{j}, [2,3,1]), nj, []);
                    else
                        % collapse the first dimension = 1
                        tmp = reshape(obj.cores{j}, nj, []);
                    end
                    T   = reshape(permute(reshape( eval(obj.oneds{j},tmp,x(xind(i),:)), nx, rj, []), [2,1,3]), rj*nx, []);
                    % how to speed up this part?
                    ii  = reshape(1:rj*nx, [], 1);
                    jj  = reshape(repmat(1:nx, rj, 1), [], 1);
                    B   = sparse(ii, jj, fx(:), rj*nx, nx);
                    %
                    fx  = T'*B;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [d,rs,ns] = size(obj)
            % Find the ranks of the TT cores and the degrees of freedom of
            % the approximation basis for each coordinate, r(0) = 1 is not
            % included.
            %
            d   = length(obj.cores);
            rs  = ones(1,d);
            ns  = ones(1,d);
            for k = 1:d-1
                rs(k) = size(obj.cores{k+1}, 1);
                ns(k) = obj.oneds{k}.num_nodes;
            end
            ns(d) = obj.oneds{d}.num_nodes;
        end
        
    end
end