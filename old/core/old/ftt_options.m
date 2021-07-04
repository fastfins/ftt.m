function options = ftt_options(varargin)
%
%Set the options for ftt constructors. If no parameters passed in, it
%returns the default options. If pass in a given option, it will modify the
%option according to the modification parameters passed in.
%
%Parameters are:
%
%method:
%  Factorisation method. Options are 'Random' and 'AMEN'. Default is 'AMEN'
%
%ng_flag:
%  A flag deciding if construct ftt of the squared root of the input function. 
%  This is useful for approximating non-negative functions. Default is false
%
%max_als:
%  Max number of ALS iterations. Default is 20
%
%err_tol:
%  Error tolerance for terminating ALS iterations. Default is 1E-4
%
%init_rank:
%  Rank of the initial tensor train, default is 10
%
%kick_rank:
%  Rank of the enrichment sample size, default is 5
%
%max_rank:
%  Max rank of each core, default is 30
%
%loc_err_tol:
%  The truncation tolerance for local SVD, default is 1E-10. The SVD is
%  truncated at singular values that is about 1E-10 relative to the largest
%  singular value
%
%int_method:
%  Interpolation methods for choosing cross indices. Default is 'MaxVol',
%  we also have 'QDEIM' and 'WEDIM', over sampling is not implemented
%
%oned_ref:
%  One d reference polynomial. See 'help setup_oned' for details. Default
%  is [], which means building the Chebyshev 2nd polynomial as oned_ref
%
%default_order:
%  Default order of the one d reference polynomial. Default is 9.
%
%Tiangang Cui, August, 2019

defaultOptions  = struct([]);
defaultMethod   = 'AMEN';
expectedMethod  = {'Random','AMEN'};
defaultNGFlag   = false;
defaultMaxALS   = 10;
defaultErrTol   = 1E-4;
defaultInitRank = 10;
defaultKickRank = 2;
defaultMaxRank  = 30;
defaultLocErrTol  = 1E-10;
defaultIntM     = 'MaxVol';
expectedIntM    = {'QDEIM','WDEIM','MaxVol'};
defaultOned     = [];
defaultOrder    = 9;
%defaultSampleSet  = [];
%defaultDebugSet = [];

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1);
%
addOptional (p,'options',  defaultOptions, @(x) isstruct(x));
addParameter(p,'method',   defaultMethod,  @(x) any(validatestring(x,expectedMethod)));
addParameter(p,'ng_flag',  defaultNGFlag,  @(x) islogical(x) && isscalar(x));
addParameter(p,'max_als',  defaultMaxALS,  validScalarPosNum);
addParameter(p,'err_tol',  defaultErrTol,  validErrTol);
addParameter(p,'init_rank',defaultInitRank,validScalarPosNum);
addParameter(p,'kick_rank',defaultKickRank,validScalarPosNum);
addParameter(p,'max_rank', defaultMaxRank, validScalarPosNum);
addParameter(p,'loc_err_tol', defaultLocErrTol,validErrTol);
addParameter(p,'int_method',  defaultIntM,    @(x) any(validatestring(x,expectedIntM)));
addParameter(p,'oned_ref', defaultOned);
addParameter(p,'defaul_order',defaultOrder,validScalarPosNum);
%addParameter(p,'sample_x', defaultSampleSet);
%addParameter(p,'debug_x',  defaultDebugSet);
%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    if ~ismember('method',cellstr(p.UsingDefaults))
        disp(['Warning: reset the method from ' options.method ' to ' tmp.method])
        options.method  = tmp.method;
    end
    if ~ismember('ng_flag',cellstr(p.UsingDefaults))
        disp(['Warning: reset the ng_flag from ' options.ng_flag ' to ' tmp.ng_flag])
        options.ng_flag = tmp.ng_flag;
    end
    if ~ismember('max_als',cellstr(p.UsingDefaults))
        options.max_als = tmp.max_als;
    end
    if ~ismember('err_tol',cellstr(p.UsingDefaults))
        options.err_tol = tmp.err_tol;
    end
    if ~ismember('init_rank',cellstr(p.UsingDefaults))
        options.init_rank = tmp.init_rank;
    end
    if ~ismember('kick_rank',cellstr(p.UsingDefaults))
        options.kick_rank = tmp.kick_rank;
    end
    if ~ismember('max_rank',cellstr(p.UsingDefaults))
        options.max_rank  = tmp.max_rank;
    end
    if ~ismember('loc_err_tol',cellstr(p.UsingDefaults))
        options.loc_err_tol = tmp.loc_err_tol;
    end
    if ~ismember('int_method',cellstr(p.UsingDefaults))
        options.int_method  = tmp.int_method;
    end
    if ~ismember('oned_ref',cellstr(p.UsingDefaults))
        options.oned_ref  = tmp.oned_ref;
    end
    if ~ismember('defaul_order',cellstr(p.UsingDefaults))
        options.defaul_order = tmp.defaul_order;
    end
end

if isempty(options.oned_ref)
    options.oned_ref = setup_oned(options.defaul_order);
end

end