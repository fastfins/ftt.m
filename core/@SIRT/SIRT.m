classdef SIRT < FTT
    % SIRT class
    %
    % SIRT Properties:
    %   marginal_direction - The direction for marginalising the FTT.
    %                 >0: marginalise from x_d to x_1
    %                 <0: marginalise from x_1 to x_d
    %   z           - Normalising constant
    %   ys          - A cell array holding squared FTT cores.
    %   ms          - A cell array for computing the marginal density.
    %   oned_cdfs   - Data structure containing information for building
    %                 one dimensional CDFs.
    %
    % FTT Methods:
    %   marginalise     - Marginalise the squared FTT.
    %   eval_irt        - X = R^{-1}(Z), where X is target random variable, 
    %                     R is Rosenblatt transport, and Z is uniform.
    %                     This function supports marginal random variables.
    %   eval_rt         - Z = R(X), where Z is uniform and X is target.
    %                     This function supports marginal random variables.
    %   eval_rt_jac     - Evaluate the Jacobian of Z = R(X).
    %   eval_cond_irt   - X_1 | X_2  = R^{-1}(Z_1;X_2), where X_2 is the
    %                     variable the Rosenblatt transport conditioned on.
    %   eval_pdf        - Evaluate the (marginal) pdf.
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example 1: (vector function outputs, m = 2):
    %
    % % Step 1: speficy the target function 
    %   func = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
    %   d = 10; % dimensionality of the input
    %
    % % Step 2: setup the Legendre basis polynomial with order 20 in 
    % % domain [0,1]
    %   poly = Legendre(20, [0,1]);
    %
    % % Step 3: use alternating energy enrichment (AMEN), default option
    %   opt = FTToption('max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-10, ...
    %           'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
    % % Optional: debug samples
    %   debug_size = 1E4;
    %   debug_x = zeros(d, debug_size);
    %   for k = 1:d
    %       debug_x(k,:) = sample_domain(poly, debug_size);
    %   end
    %
    % % Step 4: build FTT
    %   tt =  FTT(func, d, poly, opt, 'debug_x', debug_x);
    %
    % % Step 5: evaluate the function and its factorisation
    %   exact   = func(debug_x);
    %   appr_tt = eval(tt, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %
    % % Step 6: round the FTT by truncation local SVD. Here 1E-4 is the 
    % % truncation threshold of each SVD (relative to the largest singular
    % % value)
    %   ttr = round(tt, 1E-4); 
    %   appr_tt = eval(ttr, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also ONED, FTTOPTION, and FTT
    
    properties
        marginal_direction
        z
        ys
        ms
        oned_cdfs
    end
    
    methods (Static)
        function T = eval_oned_core_213(oned, core, x)
            rkm = size(core, 1);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval(oned, reshape(permute(core, [2,1,3]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_213_deri(oned, core, x)
            rkm = size(core, 1);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval_deri(oned, reshape(permute(core, [2,1,3]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_231(oned, core, x)
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval(oned, reshape(permute(core, [2,3,1]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        function T = eval_oned_core_231_deri(oned, core, x)
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = evalderi(oned, reshape(permute(core, [2,3,1]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        % [z,ys,ms] = cumint(ftt, dir)
        % Marginalise the pdf represented by ftt dimension by dimension
        
    end
    
    methods
        z = eval_rt(obj, r)
        % Evaluate squared RT z = T(r), where z is uniform and r is target r.v.

        [r,f] = eval_irt(obj, z)
        % Evaluate squared IRT r = T(z), where z is uniform

        [r,f] = eval_cond_irt(obj, x, z)
        % Using SIRT to draw conditional samples from the target pdf 
        % approximated by FTT

        J = eval_rt_jac(firt, r, z)
        % Evaluate the jacobian of the squared RT z = T(r), where z is 
        % uniform and r is target r.v.

        fx = eval_pdf(obj, x)
        % Evaluate the marginalise pdf represented by ftt

        obj = marginalise(obj, dir) 
        % Marginalise the pdf represented by ftt dimension by dimension
        
        function obj = SIRT(func, d, arg, varargin)
            % Setup data structure used for IRT
            
            defaultErrTol = 1E-8;
            %
            obj@FTT(func, d, arg, varargin{:})
            %
            if strcmp(obj.opt.tt_method, 'amen')
                obj = round(obj);
            end
            obj.oned_cdfs = cell(size(obj.oneds));
            for i = 1:length(obj.oneds)
                obj.oned_cdfs{i} = CDFconstructor(obj.oneds{i}, defaultErrTol);
            end
        end
    end
    
end