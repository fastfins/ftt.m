function [ftt, options] = build_ftt(func, d, ftt, varargin)
%Construct function tensor train for a function mapping from R^d to R^m.
%
%%%%%
%Inputs: 
%
%func:
%  A function (R^d to R^m) that take inputs as a d x n matrix and returns 
%  m x n vector
%
%d: 
%  dimension of the input variable
%
%ftt: 
%  An existing ftt, used as the initial guess
%
%options:
%  Contains options for building the ftt: 
%
%  method:
%    Factorisation method. Options are 'Random' and 'AMEN'. Default is 'AMEN'
%
%  ng_flag:
%    A flag deciding if construct ftt of the squared root of the input function. 
%    This is useful for approximating non-negative functions. Default is false
%
%  max_als:
%    Max number of ALS iterations. Default is 20
%
%  err_tol:
%    Error tolerance for terminating ALS iterations. Default is 1E-4
%
%  init_rank:
%    Rank of the initial tensor train, default is 10
%
%  kick_rank:
%    Rank of the enrichment sample size, default is 5
%
%  max_rank:
%    Max rank of each core, default is 30
%
%  loc_err_tol:
%    The truncation tolerance for local SVD, default is 1E-10. The SVD is
%    truncated at singular values that is about 1E-10 relative to the 
%    largest singular value
%
%  deim:
%    DEIM methods for choosing cross indices. Default is the QR. 'QDEIM'
%
%  oned_ref:
%    One d reference polynomial. See 'help setup_oned' for details. Default
%    is [], which means building the Chebyshev 2nd polynomial as oned_ref
%
%  default_order:
%    Default order of the one d reference polynomial. Default is 9.
%
%%%%%
%Other parameters:
%sample_x:
%  A sample set used for building the TT. It needs to be 
%  d x (init_rank + kick_rank*max_als)
%
%debug_x
%  A sample set used for debugging the TT by comparing with actual function 
%  evaluation. It needs to be d x n
%
%%%%%
%Outputs:
%
%ftt:
%  the function tensor train, it contains:
%
%  direction:   
%    ALS direction, >0: built from left to right
%                   <0: built from right to left
%
%  oneds{k} data structure:
%    contains all the information for building one d polynomial representation
%
%  cores{k}:
%    noal values of one d functions the dimension of the current core is 
%    organised as previous rank x oned{k}.num_nodes x current rank
%
%  interp_x{x}:     interpolation coordinates
%
%  ng_flag:         if the squared root of the function is factorised
%
%%%%%
%Example 1 (vector function outputs, m = 2):
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
%%%%%
%Example 2:
%
% Type 'help eval_irt' to check how to use ftt to generate random
% variables from a probability density function
%
%Tiangang Cui, August, 2019

defaultOptions   = struct([]);
defaultSampleSet = [];
defaultDebugSet  = [];

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
addRequired(p,'func',@(x) isa(x, 'function_handle'));
%
addRequired(p,'d',validScalarPosNum);
%
addRequired(p,'ftt');
%
addOptional(p,'options',defaultOptions, @(x) isstruct(x));
%
addParameter(p,'sample_x', defaultSampleSet);
addParameter(p,'debug_x',  defaultDebugSet);
%
p.KeepUnmatched = false;
parse(p,func,d,ftt,varargin{:});

options     = p.Results.options;
sample_x    = p.Results.sample_x;
debug_x     = p.Results.debug_x;


%at here, we need to get func, d, m, ftt, options, sample_x, and debug_x
if isempty(options)
    options = ftt_options();
end

% if the ftt is already constructed, flip the direction to renew
% start from new otherwise
if ~isempty(ftt)
    ftt.direction = -ftt.direction; 
else
    % copy properties
    ftt.direction = 1;
    
    % one D interpolation basis
    ftt.oneds = cell(d,1);
    if iscell(options.oned_ref)
        if length(options.oned_ref) ~= d
            error('Error: number of input reference polynomials does not match dimension')
        end
        for k = 1:d
            ftt.oneds{k} = options.oned_ref{k};
        end
    else
        for k = 1:d
            ftt.oneds{k} = options.oned_ref;
        end
    end

    % nested interpolation points
    ftt.cores = cell(d,1);
    
    % use the user provided sample set 
    if size(sample_x, 2) < options.init_rank
        disp('Not enough number of user provided samples to initialise ftt')
        sample_x = [];
        ftt.interp_x{d} = sample_oned_domain(ftt.oneds{d}, options.init_rank);
        for k = (d-1):-1:1
            ftt.interp_x{k} = [ftt.interp_x{k+1}; sample_oned_domain(ftt.oneds{k}, options.init_rank)];
        end
    else
        for k = d:-1:1
            ftt.interp_x{k} = sample_x(k:d, 1:options.init_rank);
        end
        sample_x(:,1:options.init_rank) = [];
    end
    %
    % determine m
    y = func(ftt.interp_x{k}(:,1));
    m = size(y,1);
    % interpolation weights
    ftt.cores{1} = zeros(1, ftt.oneds{1}.num_nodes, options.init_rank, m);
    for k = 2:d-1
        ftt.cores{k} = zeros(options.init_rank, ftt.oneds{k}.num_nodes, options.init_rank);
    end
    ftt.cores{d} = zeros(options.init_rank, ftt.oneds{d}.num_nodes, 1);
    
    % initialise AMEN
    if strcmp(options.method, 'AMEN')
        if size(sample_x, 2) < options.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            sample_x = [];
            ftt.res_x{d} = sample_oned_domain(ftt.oneds{d}, options.kick_rank);
            for k = (d-1):-1:1
                ftt.res_x{k} = [ftt.res_x{k+1}; sample_oned_domain(ftt.oneds{k}, options.kick_rank)];
            end
        else
            % nested interpolation points for res
            for k = d:-1:1
                ftt.res_x{k} = sample_x(k:d, 1:options.kick_rank);
            end
            sample_x(:,1:options.kick_rank) = [];
        end
        % res weights, used on the right
        for k = 1:d
            ftt.res_w{k} = randn(options.init_rank, options.kick_rank);
        end
    end
end

ftt.ng_flag = options.ng_flag;

switch options.method
    case{'Random'}
        ftt = build_ftt_rand(func, ftt, options, sample_x, debug_x);
    case{'AMEN'}
        ftt = build_ftt_amen(func, ftt, options, debug_x);
end

end

