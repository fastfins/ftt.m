classdef FTToption
    % Set the options for constructing functional tensor train
    %
    % FTToption Properties:
    %
    %   sqrt_flag   - a flag indicating if compute the sqrt of the input 
    %                 function. Default: false
    %
    %   max_rank    - max rank of each core, default is 20
    %
    %   init_rank   - rank of the initial tensor train, default is 10
    %
    %   kick_rank   - rank of the enrichment sample size, default is 3
    %
    %   max_als     - max number of ALS iterations. Default is 5
    %
    %   err_tol     - tolerance for terminating ALS. Default is 1E-2
    %
    %   loc_err_tol - truncation tolerance of local SVD, default is 1E-10.
    %                 The SVD is truncated at singular values that is about
    %                 1E-10 relative to the largest singular value
    %
    %   method      - construction method. Default is option is 'amen'. 
    %                 Options are 'amen' and 'random'.
    %
    %   int_method  - interpolation method for choosing cross indices. 
    %                 Default is 'MaxVol', over sampling is not implemented
    %
    % FTToption Methods:
    %
    %    FTToption  - constructor. If no parameter is passed in, it returns 
    %                 the default obj.
    %
    %    update     - update options of an object
    %
    
    properties
        sqrt_flag
        loc_err_tol
        max_rank
        init_rank
        kick_rank
        max_als
        err_tol
        method
        int_method
    end
    
    properties (Access = private, Constant = true)
        defaultSQRTFlag = false;
        defaultMaxALS   = 5;
        defaultErrTol   = 1E-2;
        defaultInitRank = 10;
        defaultKickRank = 3;
        defaultMaxRank  = 20;
        defaultLocErrTol  = 1E-10;
        defaultMethod   = 'amen';
        expectedMethod  = {'random','amen'};
        defaultIntM     = 'MaxVol';
        expectedIntM    = {'QDEIM','WDEIM','MaxVol'};
    end
    
    methods
        function obj = FTToption(varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1);
            %
            addParameter(p,'sqrt_flag',obj.defaultSQRTFlag,@(x) islogical(x) && isscalar(x));
            addParameter(p,'max_als',  obj.defaultMaxALS,  validScalarPosNum);
            addParameter(p,'err_tol',  obj.defaultErrTol,  validErrTol);
            addParameter(p,'init_rank',obj.defaultInitRank,validScalarPosNum);
            addParameter(p,'kick_rank',obj.defaultKickRank,validScalarPosNum);
            addParameter(p,'max_rank', obj.defaultMaxRank, validScalarPosNum);
            addParameter(p,'loc_err_tol', obj.defaultLocErrTol,validErrTol);
            addParameter(p,'method',  obj.defaultMethod, ...   
                @(x) any(validatestring(x,obj.expectedMethod)));
            addParameter(p,'int_method',  obj.defaultIntM, ...   
                @(x) any(validatestring(x,obj.expectedIntM)));
            %
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            tmp = p.Results;
            %
            obj.sqrt_flag = tmp.sqrt_flag;
            obj.max_als   = tmp.max_als;
            obj.err_tol   = tmp.err_tol;
            obj.init_rank = tmp.init_rank;
            obj.kick_rank = tmp.kick_rank;
            obj.max_rank  = tmp.max_rank;
            obj.loc_err_tol = tmp.loc_err_tol;
            obj.method      = tmp.method;
            obj.int_method  = tmp.int_method;
        end
        
        function obj = update(obj, varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1);
            %
            addParameter(p,'sqrt_flag',obj.defaultSQRTFlag,@(x) islogical(x) && isscalar(x));
            addParameter(p,'max_als',  obj.defaultMaxALS,  validScalarPosNum);
            addParameter(p,'err_tol',  obj.defaultErrTol,  validErrTol);
            addParameter(p,'init_rank',obj.defaultInitRank,validScalarPosNum);
            addParameter(p,'kick_rank',obj.defaultKickRank,validScalarPosNum);
            addParameter(p,'max_rank', obj.defaultMaxRank, validScalarPosNum);
            addParameter(p,'loc_err_tol', obj.defaultLocErrTol,validErrTol);
            addParameter(p,'method',  obj.defaultMethod, ...   
                @(x) any(validatestring(x,obj.expectedMethod)));
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
            if ~ismember('err_tol',cellstr(p.UsingDefaults))
                obj.err_tol = tmp.err_tol;
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
            if ~ismember('loc_err_tol',cellstr(p.UsingDefaults))
                obj.loc_err_tol = tmp.loc_err_tol;
            end
            if ~ismember('method',cellstr(p.UsingDefaults))
                obj.method = tmp.method;
            end
            if ~ismember('int_method',cellstr(p.UsingDefaults))
                obj.int_method = tmp.int_method;
            end
        end
    end
end