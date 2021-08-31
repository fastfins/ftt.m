classdef DIRT
    % DIRT class
    %
    % DIRT Properties:    
    %   diag        - Type of the diagonal transformation in DIRT
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
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example:
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also SIRT
    
    properties
        diag
        logz
        irts
        method
        %
        d
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
        [r,f] = eval_irt(obj, z)
        % Evaluate DIRT r = T(z), where z follows some general reference
        
        [r,f] = eval_cirt(obj, x, z)
        % Evaluate the conditional DIRT
        
        z = eval_rt(obj, r)
        % Evaluate deep RT z = T(r), where z is reference and r is target r.v.
        
        obj = build(obj, func)
        % building DIRT using given temperatures
                
        f = ratio_fun(obj, func, beta_p, beta, z)
        % ratio function for building DIRT
               
        function obj = DIRT(func, d, beta, varargin)
            % Call FTT constructor to build the FTT and setup data
            % structures for SIRT. Need to run marginalise after this.
            % parsing inputs
            defaultPoly     = Lagrange1(20, [0,1]);
            defaultOption   = FTToption();
            defaultMethod   = 'Aratio';
            expectedMethod  = {'Eratio','Aratio'};
            defaultDiag     = uniformMap();
            defaultMLayers  = 50;
            %defaultQMCFlag  = false;
            defaultNSamples = 1E3;
            defaultNDebugs  = 1E3;
            defaultMinBeta  = 1E-6;
            defaultESSTol   = 0.2;
            defaultBetas    = [];
            
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && all(x > 0);
            %
            addRequired(p, 'func',  @(x) isa(x, 'function_handle'));
            addRequired(p, 'd',     validScalarPosNum);
            addRequired(p, 'beta',  @(x) isnumeric(x) && (x > 0));
            addOptional(p, 'poly',  defaultPoly);
            addOptional(p, 'option',defaultOption);
            addParameter(p,'method',defaultMethod, ...
                @(x) any(validatestring(x,expectedMethod)));
            addParameter(p,'diag',  defaultDiag);
            addParameter(p,'max_layers',defaultMLayers, validScalarPosNum);
            %addParameter(p,'qmc_flag',  defaultQMCFlag, @(x) islogical(x) && isscalar(x));
            addParameter(p,'n_samples', defaultNSamples,validScalarPosNum);
            addParameter(p,'n_debugs',  defaultNDebugs, validScalarPosNum);
            addParameter(p,'min_beta',  defaultMinBeta, validScalarPosNum);
            addParameter(p,'ess_tol',   defaultESSTol,  validScalarPosNum);
            addParameter(p,'betas', defaultBetas);
            %
            p.KeepUnmatched = false;
            parse(p,func,d,beta,varargin{:});
            %
            oned = p.Results.poly;
            sirt_opt = p.Results.option;
            obj.d = d;
            obj.method = p.Results.method;
            obj.diag = p.Results.diag;
            obj.max_layers = p.Results.max_layers;
            %obj.qmc_flag = p.Results.qmc_flag;
            obj.n_samples = p.Results.n_samples;
            obj.n_debugs = p.Results.n_debugs;
            %
            obj.min_beta = p.Results.min_beta;
            obj.ess_tol = p.Results.ess_tol;
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
            obj = build(obj, func, oned, sirt_opt);
        end
    end
    
end