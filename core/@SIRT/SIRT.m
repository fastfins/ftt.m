classdef SIRT < FTT
    % SIRT class
    %
    % SIRT Properties:
    %   int_dir     - The direction for marginalising the FTT.
    %                 >0: marginalise from x_d to x_1
    %                 <0: marginalise from x_1 to x_d
    %   z           - Normalising constant
    %   ys          - A cell array holding squared FTT cores.
    %   ms          - A cell array for computing the marginal density.
    %   oned_cdfs   - One dimensional bases for building CDFs.
    %
    % FTT Methods:
    %   marginalise - Marginalise the squared FTT.
    %   eval_pdf    - Evaluate the normalised (marginal) pdf.
    %   eval_irt    - X = R^{-1}(Z), where X is the target random variable, 
    %                 R is the Rosenblatt transport, and Z is the uniform
    %                 random variable. 
    %               * This function can map marginal random variables.
    %   eval_cirt   - Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly 
    %                 follow the target represented by SIRT, Z is uniform. 
    %               * This function cannot handle marginal random variables.
    %   eval_rt     - Z = R(X), where Z is uniform and X is target.
    %               * This function can map marginal random variables.
    %   eval_rt_jac - Evaluate the Jacobian of Z = R(X).
    %               * This function cannot handle marginal random variables.
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example: 
    %
    % % Setup a d-dimensional Gausssian density with correlation controlled
    % % by a.
    %   d = 20; a = 0.5;
    %   A = diag(-sqrt(1-a^2)*ones(d-1,1), -1) + eye(d);
    %   D = diag([1, a*ones(1,d-1)]);
    %   B = D\A;
    %   Q = B'*B;   % precision matrix
    %   C = inv(Q); % covariance matrix
    %   z = sqrt((2*pi)^d/det(Q)); % normalising constant
    % % The joint distribution, unnormalised
    %   joint = @(x) exp(-0.5*sum((B*x).^2,1));
    %
    % % Build the SIRT using Fourier basis
    %   pol = Fourier(20, [-5,5]);
    %   opt = FTToption('als_tol',1E-4,'max_rank',20,'sqrt_flag',true);
    %   dx  = B\randn(d, 1E4);
    %   sx  = B\randn(d, 1E3);
    %   irt =  SIRT(joint, d, pol, opt, 'debug_x', dx, 'sample_x', sx);
    %
    %%%%%%%%%%%%%%%%%
    %
    % % Task 1: inverse Rosenblatt and Rosenblatt transports
    %   Z   = rand(d, 1E4);     % uniform random seeds
    %   [X,f] = eval_irt(irt,Z);% run IRT, f is the joint density
    %   U   = eval_rt (irt,X);  % U should be the same as Z
    %   p   = eval_pdf(irt,X);  % should be the same as f
    %   fe  = joint(X);         % exact pdf
    %   norm(Z - U, 'fro')
    %   norm(f - p, 'fro')
    %   figure; plot(fe, f/z, '.')
    %   figure; plot(C - cov(X')); % error in the estimated covariance
    %
    %%%%%%%%%%%%%%%%%
    %
    % % Task 2: inverse Rosenblatt and Rosenblatt transports for marginal
    % % random variables.
    % % The marginal distribution, unnormalised
    %   marginal = @(x, ind) exp(-0.5*sum((C(ind,ind)\x).*x,1));
    % % normalising constant of the marginal
    %   marginal_z = @(ind) sqrt(det(C(ind,ind))*(2*pi)^length(ind));
    %
    % % Case 2.1: Marginal for X_1, ..., X_k 
    % %     - need to run the marginal function with int_dir=1 (default)
    %   if irt.int_dir ~= 1, irt = marginalise(irt, 1); end
    %   ind = 1:8; % the leftover coordinates after marginalisation
    %   [Xm,f] = eval_irt(irt, Z(ind,:));
    %   p = eval_pdf(irt, Xm);
    %   Um = eval_rt(irt, Xm);
    %   norm(Z(ind,:) - Um, 'fro')
    %   norm(f - p, 'fro')
    %   fe = marginal(Xm, ind)/marginal_z(ind);
    %   figure; plot(fe , f, '.');
    %   figure; plot(C(ind, ind) - cov(Xm'))
    %
    % % Case 2.2: Marginal for X_{k+1}, ..., X_d 
    % %     - need to run the marginal function with int_dir=-1
    %   if irt.int_dir ~= -1, irt = marginalise(irt, -1); end
    %   ind = 13:d; % the leftover coordinates after marginalisation
    %   [Xm,f] = eval_irt(irt, Z(ind,:));
    %   p = eval_pdf(irt, Xm);
    %   Um = eval_rt(irt, Xm);
    %   norm(Z(ind,:) - Um, 'fro')
    %   norm(f - p, 'fro')
    %   fe = marginal(Xm, ind)/marginal_z(ind);
    %   figure; plot(fe , f, '.');
    %   figure; plot(C(ind, ind) - cov(Xm'))
    %
    %%%%%%%%%%%%%%%%%
    %
    % % Task 3: inverse Rosenblatt transport for conditional random variables.
    % % Case 3.1: X_{>=k} | x_{<k} 
    % %     - need to run the marginal function with int_dir=1 (default)
    %   indx = 1:8; % index of variables that will be conditioned on
    %   indy = 9:d; % index of conditional variables 
    %   X = B\randn(d,1);
    %   X = X(indx);
    %   % conditional mean
    %   my = C(indy,indx)*(C(indx,indx)\X); 
    %   % conditional covariance
    %   Cy = C(indy,indy) - C(indy,indx)*(C(indx,indx)\C(indx,indy));
    %   Zy = Z(indy,:);
    %   if irt.int_dir ~= 1, irt = marginalise(irt, 1); end
    %   [Y,f] = eval_cirt(irt, X, Zy);
    %   fe = joint([repmat(X,1,size(Zy,2)); Y])/marginal(X,indx)*marginal_z(indx)/z;
    %   figure; plot(fe , f, '.');
    %   figure; plot(Cy - cov(Y'))
    %   figure; plot(my - mean(Y,2))
    %
    % % Case 3.2: X_{<=k} | x_{>k}
    % %     - need to run the marginal function with int_dir=-1
    %   indx = 9:d; % index of variables that will be conditioned on
    %   indy = 1:8; % index of conditional variables 
    %   X = B\randn(d,1);
    %   X = X(indx);
    %   % conditional mean
    %   my = C(indy,indx)*(C(indx,indx)\X); 
    %   % conditional covariance
    %   Cy = C(indy,indy) - C(indy,indx)*(C(indx,indx)\C(indx,indy));
    %   Zy = Z(indy,:);
    %   if irt.int_dir ~= -1, irt = marginalise(irt, -1); end
    %   [Y,f] = eval_cirt(irt, X, Zy);
    %   fe = joint([Y; repmat(X,1,size(Zy,2))])/marginal(X,indx)*marginal_z(indx)/z;
    %   figure; plot(fe , f, '.');
    %   figure; plot(Cy - cov(Y'))
    %   figure; plot(my - mean(Y,2))
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also ONED, ONEDCDF, FTTOPTION, and FTT
    
    properties
        int_dir
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
            tmp = eval_deri(oned, reshape(permute(core, [2,3,1]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        % [z,ys,ms] = cumint(ftt, dir)
        % Marginalise the pdf represented by ftt dimension by dimension
        
    end
    
    methods
        [r,f] = eval_irt(obj, z)
        % Evaluate squared IRT r = T(z), where z is uniform

        [r,f] = eval_cirt(obj, x, z)
        % Using SIRT to draw conditional samples from the target pdf 
        % approximated by FTT

        z = eval_rt(obj, r)
        % Evaluate squared RT z = T(r), where z is uniform and r is target r.v.
        
        J = eval_rt_jac(firt, r, z)
        % Evaluate the jacobian of the squared RT z = T(r), where z is 
        % uniform and r is target r.v.

        fx = eval_pdf(obj, x)
        % Evaluate the marginalise pdf represented by ftt

        obj = marginalise(obj, dir) 
        % Marginalise the pdf represented by ftt dimension by dimension
        
        function obj = SIRT(func, d, arg, varargin)
            % Call FTT constructor to build the FTT and setup data 
            % structures for SIRT. Need to run marginalise after this.
            
            %
            obj@FTT(func, d, arg, varargin{:})
            %
            if strcmp(obj.opt.tt_method, 'amen')
                obj = round(obj);
            end
            %
            obj.oned_cdfs = cell(size(obj.oneds));
            for i = 1:length(obj.oneds)
                obj.oned_cdfs{i} = CDFconstructor(obj.oneds{i}, obj.opt.cdf_tol);
            end
            obj = marginalise(obj);
        end
    end
    
end