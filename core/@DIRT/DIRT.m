classdef DIRT
    % DIRT class
    %
    % DIRT Properties:
    %   int_dir     - The direction for marginalising the FTT.
    %                 >0: marginalise from x_d to x_1
    %                 <0: marginalise from x_1 to x_d
    %   method      - Type of DIRT construction, choose either 'Eratio'
    %                 or 'Aratio', default is 'Aratio'
    %   diag        - Type of the diagonal transformation in DIRT
    %   z           - Normalising constant
    %   irts        - A cell array of IRTs
    %   betas       - A vector of temperatures used in DIRT construction
    %   n_layers    - Number of layers
    %   max_layers  - The maximum number of layers
    %   qmc_flag    - A flag to determine if use QMC samples
    %   n_samples   - 
    %
    % DIRT Methods:
    %   marginalise - Marginalise the DIRT.
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
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example: 
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also SIRT
    
    properties
        int_dir
        diag
        z
        irts
        beta tempering
    end
    
    methods (Static)

    end
    
    methods
        [r,f] = eval_irt(obj, z)
        % Evaluate squared IRT r = T(z), where z is uniform

        [r,f] = eval_cirt(obj, x, z)
        % Using SIRT to draw conditional samples from the target pdf 
        % approximated by FTT

        z = eval_rt(obj, r)
        % Evaluate squared RT z = T(r), where z is uniform and r is target r.v.
        
        % J = eval_rt_jac(firt, r, z)
        % Evaluate the jacobian of the squared RT z = T(r), where z is 
        % uniform and r is target r.v.

        fx = eval_pdf(obj, x)
        % Evaluate the marginalise pdf represented by ftt

        obj = marginalise(obj, dir) 
        % Marginalise the pdf represented by ftt dimension by dimension
        
        obj = tempering(func, d, beta, pol, opt)
        % building DIRT using given temperatures
        
        obj = adapt_tempering(func, d, beta, pol, opt)
        % building DIRT using adaptive temperatures

        f = ratio_fun(obj, func, beta_p, beta, v)
        % ratio function for building DIRT

        function obj = DIRT(func, d, beta, varargin)
            % Call FTT constructor to build the FTT and setup data 
            % structures for SIRT. Need to run marginalise after this.
            % parsing inputs
            defaultPoly     = Lagrange1(20, [0,1]);
            defaultOption   = FTToption();
            defaultMethod   = 'Aratio';
            expectedMethod  = {'Eratio','Aratio'};
            defaultDiag     = 'uniform';
            expectedDiag    = {'uniform','normal'};
            defaultMLayers  = 50;
            defaultQMCFlag  = false;
            
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
            addParameter(p,'diag',  defaultDiag, ...   
                @(x) any(validatestring(x,expectedDiag)));
            addParameter(p,'max_layers',defaultMLayers, validScalarPosNum);
            addParameter(p,'qmc_flag',  defaultQMCFlag, @(x) islogical(x) && isscalar(x));
            %
            p.KeepUnmatched = false;
            parse(p,func,d,beta,varargin{:});
            %
            pol = p.Results.poly;
            opt = p.Results.option;
            obj.method = p.Results.method;
            obj.diag = p.Results.diag;
            obj.max_layers = p.Results.max_layers;
            obj.qmc_flag = p.Results.qmc_flag;
            %
            obj = build(obj, func, d, beta, pol, opt);
        end
    end
    
end