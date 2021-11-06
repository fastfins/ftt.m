classdef FTToption
    % Set the options for constructing functional tensor train
    %
    % FTToption Properties:
    %   sqrt_flag   - A flag indicating if compute the sqrt of the input 
    %                 function. Default is false.
    %   max_rank    - Max rank of each core, default is 20.
    %   init_rank   - Rank of the initial tensor train, default is 10.
    %   kick_rank   - Rank of the enrichment sample size, default is 2.
    %   max_als     - Max number of ALS iterations. Default is 4.
    %   als_tol     - Tolerance for terminating ALS. Default is 1E-2.
    %   local_tol   - Truncation tolerance of local SVD, default is 1E-10.
    %                 The SVD is truncated at singular values that is about
    %                 1E-10 relative to the largest singular value.
    %   cdf_tol     - Tolerance for evaluating the inverse CDF function.
    %                 Used by SIRT. Default is 1E-10.
    %   tt_method   - Construction method. Default is option is 'amen'. 
    %                 Options are 'amen' and 'random'.
    %   int_method  - Interpolation method for choosing cross indices. 
    %                 Default is 'MaxVol'.
    %
    % FTToption Methods:
    %    FTToption  - Constructor. If no parameter is passed in, it returns 
    %                 the default obj.
    %    update     - Update property 
    %
    
    properties
        sqrt_flag
        local_tol
        max_rank
        init_rank
        kick_rank
        max_als
        als_tol
        cdf_tol
        tt_method
        int_method
    end
    
    properties (Access = private, Constant = true)
        defaultSQRTFlag = false;
        defaultMaxALS   = 4;
        defaultALSTol   = 1E-2;
        defaultInitRank = 10;
        defaultKickRank = 2;
        defaultMaxRank  = 20;
        defaultLocTol   = 1E-10;
        defaultCDFTol   = 1E-10;
        defaultTTMethod = 'amen';
        expectedTTMethod  = {'random','amen'};
        defaultIntM     = 'MaxVol';
        expectedIntM    = {'QDEIM','WDEIM','MaxVol'};
    end
    
    methods
        function obj = FTToption(varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            addParameter(p,'sqrt_flag',obj.defaultSQRTFlag,@(x) islogical(x) && isscalar(x));
            addParameter(p,'max_als',  obj.defaultMaxALS,  validScalarPosNum);
            addParameter(p,'als_tol',  obj.defaultALSTol,  validErrTol);
            addParameter(p,'init_rank',obj.defaultInitRank,validScalarPosNum);
            addParameter(p,'kick_rank',obj.defaultKickRank,validScalarPosNum);
            addParameter(p,'max_rank', obj.defaultMaxRank, validScalarPosNum);
            addParameter(p,'local_tol',obj.defaultLocTol,  validErrTol);
            addParameter(p,'cdf_tol',  obj.defaultCDFTol,  validErrTol);
            addParameter(p,'tt_method',obj.defaultTTMethod, ...   
                @(x) any(validatestring(x,obj.expectedTTMethod)));
            addParameter(p,'int_method',  obj.defaultIntM, ...   
                @(x) any(validatestring(x,obj.expectedIntM)));
            %
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            tmp = p.Results;
            %
            obj.sqrt_flag = tmp.sqrt_flag;
            obj.max_als   = tmp.max_als;
            obj.als_tol   = tmp.als_tol;
            obj.init_rank = tmp.init_rank;
            obj.kick_rank = tmp.kick_rank;
            obj.max_rank  = tmp.max_rank;
            obj.local_tol = tmp.local_tol;
            obj.cdf_tol   = tmp.cdf_tol;
            obj.tt_method = tmp.tt_method;
            obj.int_method  = tmp.int_method;
        end
        
        function obj = update(obj, varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            addParameter(p,'sqrt_flag',obj.defaultSQRTFlag,@(x) islogical(x) && isscalar(x));
            addParameter(p,'max_als',  obj.defaultMaxALS,  validScalarPosNum);
            addParameter(p,'als_tol',  obj.defaultALSTol,  validErrTol);
            addParameter(p,'init_rank',obj.defaultInitRank,validScalarPosNum);
            addParameter(p,'kick_rank',obj.defaultKickRank,validScalarPosNum);
            addParameter(p,'max_rank', obj.defaultMaxRank, validScalarPosNum);
            addParameter(p,'local_tol',obj.defaultLocTol,  validErrTol);
            addParameter(p,'tt_method',obj.defaultTTMethod, ...   
                @(x) any(validatestring(x,obj.expectedTTMethod)));
            addParameter(p,'int_method',  obj.defaultIntM,... 
                @(x) any(validatestring(x,obj.expectedIntM)));
            %
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            tmp = p.Results;
            %
            if ~ismember('sqrt_flag',cellstr(p.UsingDefaults))
                disp(['Warning: reset the sqrt_flag from ' obj.sqrt_flag ' to ' tmp.sqrt_flag])
                obj.sqrt_flag = tmp.sqrt_flag;
            end
            if ~ismember('max_als',cellstr(p.UsingDefaults))
                obj.max_als = tmp.max_als;
            end
            if ~ismember('als_tol',cellstr(p.UsingDefaults))
                obj.als_tol = tmp.als_tol;
            end
            if ~ismember('init_rank',cellstr(p.UsingDefaults))
                obj.init_rank = tmp.init_rank;
            end
            if ~ismember('kick_rank',cellstr(p.UsingDefaults))
                obj.kick_rank = tmp.kick_rank;
            end
            if ~ismember('max_rank',cellstr(p.UsingDefaults))
                obj.max_rank = tmp.max_rank;
            end
            if ~ismember('local_tol',cellstr(p.UsingDefaults))
                obj.local_tol = tmp.local_tol;
            end
            if ~ismember('cdf_tol',cellstr(p.UsingDefaults))
                obj.cdf_tol = tmp.cdf_tol;
            end
            if ~ismember('tt_method',cellstr(p.UsingDefaults))
                obj.tt_method = tmp.tt_method;
            end
            if ~ismember('int_method',cellstr(p.UsingDefaults))
                obj.int_method = tmp.int_method;
            end
        end
    end
end