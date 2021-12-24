classdef DIRT
    % DIRT class
    %
    % DIRT Properties:    
    %   ref         - The diagonal transformation and reference in DIRT
    %   logz        - Log normalising constant
    %   irts        - A cell array of IRTs
    %   method      - Type of DIRT construction, choose either 'Eratio'
    %                 or 'Aratio', default is 'Aratio'
    %   d           - dimension of the inpuit parameters
    %   n_layers    - Number of layers
    %   max_layers  - The maximum number of layers
    %   n_samples   - Number of samples used in the temperature adaptation 
    %               - and in the FTT construction
    %   n_debugs    - Number of debug samples
    %   adapt_beta  - A flag indicating if adaptive temperature is used
    %   min_beta    - Minimum temperature
    %   ess_tol     - Tolerance for increasing the temperature
    %   betas       - List of temperatures
    %   n_evals     - Number of function evaluations in TT cross.
    %
    % DIRT Methods:
    %   eval_irt    - X = R^{-1}(Z), where X is the target random variable,
    %                 R is the Rosenblatt transport, and Z is the uniform
    %                 random variable.
    %               * This function can map marginal random variables.
    %   eval_cirt   - Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly
    %                 follow the target represented by SIRT, Z is uniform.
    %               * This function cannot handle marginal random variables.
    %   eval_rt     - Z = R(X), where Z is uniform and X is target.
    %               * This function can map marginal random variables.
    %   backtrack   - Compute the gradient of f(T(z)) w.r.t. z
    %   ratio_fun   - The ratio function used in DIRT construction
    %   random      - Generate random variables
    %   sobol       - Generate transformed sobol points
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example: 
    %
    % % Use the following function in the example
    %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   function [mllkd, mlp] = fun_banana(u, data, sigma)
    %   F = log((1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2);
    %   mllkd = sum((F-data).^2,1)/(2*sigma^2);
    %   mlp = 0.5*sum(u.^2,1);
    %   end
    %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Parameters for the target function:
    %   data=[3;5]; sigma=0.3;
    %
    % % Use the default option, Gaussian reference with 2nd order Langrange 
    % % basis, only need to provide the density funtion, dimension of
    % % random variables and domain
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,[-4,4]);
    %
    % % Case 1: using a uniform reference measure with Lagrange basis
    % %
    % % Define the map to the reference measure
    %   ref = UniformReference();
    % % Define two basis functions, the first one is for discretising the 
    % % 0th level bridging density and the second one is for discretising 
    % % w.r.t. the reference measure.
    %   poly = {Lagrange1(60,[-4,4]), Lagrange1(60,ref.domain)};
    % % Define options for running ALS
    %   opt = FTToption('max_als',4,'init_rank',40,'max_rank', 50);
    % % Build DIRT using the approximate ratio method. The second argument
    % % is the dimenion of the parameters, 'ess_tol' is used for adaptation
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,poly,ref,opt,...
    %       'min_beta',1E-3,'ess_tol',0.8,'method','Aratio');
    %
    % % Case 2: using a uniform reference measure with Legendre basis
    % %
    % % Define the map to the reference measure
    %   ref = UniformReference();
    % % Define two basis functions, the first one is for discretising the 
    % % 0th level bridging density and the second one is for discretising 
    % % w.r.t. the reference measure.
    %   poly = {Legendre(60, [-4, 4]), Legendre(40, ref.domain)};
    % % Define options for running ALS
    %   opt = FTToption('max_als',4,'init_rank',40,'max_rank', 50);
    % % Build DIRT using the eaxct ratio method. Here the temperatures are
    % % prespecifed via the parameter 'betas'.
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,poly,ref,opt,...
    %       'betas',2.^(-9:0),'method','Eratio');
    %
    % % Case 3: using a Gaussian reference measure with Fourier
    % %
    % % Define the map to the reference measure
    %   ref = GaussReference(0,1,[-5,5]);
    % % Define two basis functions, the first one is for discretising the 
    % % 0th level bridging density and the second one is for discretising 
    % % w.r.t. the reference measure.
    %   poly = {Fourier(30,[-5,5]), Fourier(30,ref.domain)};
    % % Define options for running ALS
    %   opt = FTToption('max_als',4,'init_rank',40,'max_rank', 50);
    % % Build DIRT using the approximate ratio method. The second argument
    % % is the dimenion of the parameters, 'ess_tol' is used for adaptation
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,poly,ref,opt,...
    %       'min_beta',1E-3,'ess_tol',0.8,'method','Aratio');
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also SIRT
    
    properties
        logz
        irts
        method
        %
        d
        ref
        %
        n_layers
        max_layers
        %
        n_samples
        n_debugs
        %
        adapt_beta
        min_beta
        ess_tol
        betas
        %
        n_evals
    end
    
    methods (Static)
        gz = backtrack(Juz, Jux, gx)
        % Evaluate the gradient of f(T(z)), where T is a DIRT
    end
    
    methods
        [r,f,gz,Juz,Jux] = eval_irt(obj, z, k)
        % Evaluate DIRT r = T(z), where z follows some general reference
        
        [r,f] = eval_cirt(obj, x, z)
        % Evaluate the conditional DIRT
        
        [z,logz] = eval_rt(obj, r, k)
        % Evaluate deep RT z = T(r), where z is reference and r is target r.v.
        
        logz = log_pdf(obj, r, k)
        
        obj = build(obj, func, oneds, sirt_opt)
        % building DIRT using given temperatures
                
        f = ratio_fun(obj, func, z, sqrt_flag)
        % ratio function for building DIRT
               
        [f,g] = pullback(obj, func, z)
        % pullback of the target density function
        
        function r = random(obj, n)
            % pseudo random samples
            z = random(obj.ref, obj.d, n);
            r = eval_irt(obj, z);
        end
        
        function r = sobol(obj, n)
            % QMC samples using Sobol sequence
            z = sobol(obj.ref, obj.d, n);
            r = eval_irt(obj, z);
        end
        
        function obj = DIRT(func, d, arg, varargin)
            % Call FTT constructor to build the FTT and setup data
            % structures for SIRT. Need to run marginalise after this.
            % parsing inputs
            defaultOption   = FTToption('max_als', 2);
            defaultRef = GaussReference();
            defaultMethod   = 'Aratio';
            expectedMethod  = {'Eratio','Aratio'};
            defaultMLayers  = 50;
            %defaultQMCFlag  = false;
            defaultNSamples = 1E3;
            defaultNDebugs  = 1E3;
            defaultMinBeta  = 1E-6;
            defaultESSTol   = 0.5;
            defaultBetas    = [];
            
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && all(x > 0);
            %
            addRequired(p, 'func',  @(x) isa(x, 'function_handle'));
            addRequired(p, 'd',     validScalarPosNum);
            %addRequired(p, 'oneds');
            addRequired(p, 'arg');
            addOptional(p, 'ref', defaultRef);
            %
            addOptional(p, 'option',defaultOption);
            addParameter(p,'method',defaultMethod, ...
                @(x) any(validatestring(x,expectedMethod)));
            addParameter(p,'max_layers',defaultMLayers, validScalarPosNum);
            %addParameter(p,'qmc_flag',  defaultQMCFlag, @(x) islogical(x) && isscalar(x));
            addParameter(p,'n_samples', defaultNSamples,validScalarPosNum);
            addParameter(p,'n_debugs',  defaultNDebugs, validScalarPosNum);
            addParameter(p,'min_beta',  defaultMinBeta, validScalarPosNum);
            addParameter(p,'ess_tol',   defaultESSTol,  validScalarPosNum);
            addParameter(p,'betas', defaultBetas);
            %
            p.KeepUnmatched = false;
            parse(p,func,d,arg,varargin{:});
            %
            sirt_opt = p.Results.option;
            obj.d = d;
            obj.ref = p.Results.ref;
            obj.method = p.Results.method;
            obj.max_layers = p.Results.max_layers;
            %obj.qmc_flag = p.Results.qmc_flag;
            obj.n_samples = p.Results.n_samples;
            obj.n_debugs = p.Results.n_debugs;
            %
            obj.min_beta = p.Results.min_beta;
            obj.ess_tol = p.Results.ess_tol;
            %
            % reconstruct oneds
            if isa(arg, 'cell')
                % This should contain at most 2 cells: for levels 0 and 1
                if (numel(arg)>2)
                    warning('poly cells 3:%d are not used and will be ignored', numel(arg));
                end
                if (isa(arg{1}, 'Oned'))
                    % Replicate the same poly, since build will assume
                    % anisotropic bases always
                    oneds{1} = repmat(arg(1), 1, d);
                elseif (isa(arg{1}, 'cell'))
                    assert(numel(arg{1})==d, 'Anisotropic level 0 poly contains %d elements instead of %d', numel(arg{1}), d);
                    oneds{1} = arg{1};
                else
                    error('wrong level 0 poly argument')
                end
                assert(isa(arg{2}, 'Oned'), 'wrong level 1 poly argument');
                oneds{2} = arg{2};
%                 for k = 1:length(arg)
%                     if ~isa(arg{k}, 'Oned')
%                         error('wrong type of argument')
%                     end
%                     oneds = arg;
%                 end
            elseif isa(arg, 'Oned')
                % iniitialize the second polynomial using the first
                oneds = cell(2,1);
                oneds{1} = repmat({arg}, 1, d);  % allow anisotropic build
                if isa(arg, 'Lagrange1')
                    oneds{2} = feval(class(arg), arg.num_elems, obj.ref.domain); 
                elseif isa(arg, 'Lagrangep')
                    oneds{2} = feval(class(arg), arg.order, arg.num_elems, obj.ref.domain);
                else
                    oneds{2} = feval(class(arg), arg.order, obj.ref.domain);
                end
            elseif isa(arg, 'numeric') && length(arg) == 2
                % arg give the domain
                oneds = cell(2,1);
                oneds{1} = repmat({Lagrangep(2, 25, arg)}, 1, d);  % allow anisotropic build
                oneds{2} = Lagrangep(2, 25, obj.ref.domain);
            else
                error('wrong type of argument')
            end
            %
            betas = p.Results.betas;
            if isempty(betas)
                obj.adapt_beta = true;
                obj.betas = zeros(obj.max_layers, 1);
            else
                obj.adapt_beta = false;
                obj.betas = sort(betas, 'ascend');
                obj.max_layers = max(obj.max_layers, length(obj.betas));
            end
            %
            obj = build(obj, func, oneds, sirt_opt);
        end
    end
    
end