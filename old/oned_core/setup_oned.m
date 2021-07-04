function [def, options] = setup_oned(order, varargin)
%Define the one dimensional polynomial basis for the tensor train.
%
%Currently support piecewise Lagrange polynomials of arbitrary orders, 
%orthogonal polynomials including Chebyshev, Jacobi, Legendre, Fourier, 
%Laguerre, and Hermite. 
% 
%It can also be extend to Sinc functions, ATM the Laguerre, Hermite and
%sinc functions may not albe to support the computation of the CDF
%functions in a numerical stable way (the indefinite integral of the squared 
%approximation is not easy to compute).
%
%%%%%
%Input: 
%
%order:    
%  Max order of the global polynomials or the max order of the local
%  Lagrange polynomial
%
%%%%%
%Optional input parameters:
%
%type:     
%  Type of the basis, options are 'Lagrange', 'Chebyshev1st', 'Chebyshev2nd', 
%  'Legendre', 'Jacobi11', 'Fourier', 'Laguerre' and 'Hermite'
%  Default type: 'Chebyshev2nd'
%
%domain: 
%  The domain of the function, default is [0,1]
%
%lag_elems:
%  Number of elements if the h-p FEM basis is used, default is 1
%
%lag_grid:
%  The grid used for defining the elements if h-p FEM basis is used,
%  default is 1.
%
%lag_nodes:
%  Type of nodes used for Lagrange basis, default is 'Jacobi', can choose
%  between 'Jacobi' and 'Chebyshev'
%  
%bc:
%  Boundary condition, can choose between 'Dirichlet' (zero at bouyndary) 
%  and 'Neumann' (zero flux at boundary), default is 'Dirichlet'
%
%ghost_size: 
%  Size of the ghost element for handling boundary condition, default is 0
%
%
%%%%%
%Outputs:   
%
%def:      
%  a data structure contains neccessary for operating with basis polynomials, 
%  type in 'help setup_lagrange' and 'help setup_orthogonal' to see details.
%
%options:
%  options used for setting up the basis polynomials
%
%%%%%
%Examples:
%
%
%Tiangang Cui, August, 2019


defaultDomain   = [0, 1];
defaultType     = 'Chebyshev2nd';
expectedType    = {'Lagrange','Chebyshev1st','Chebyshev2nd','Legendre',...
    'Jacobi11','Laguerre','Hermite','Fourier'};
%
% options for Lagrange polynomials
defaultLagElems = 1;
defaultLagGrid  = [];
defaultLagNodes = 'Jacobi';
expectedLagNodes    = {'Jacobi','Chebyshev'};
defaultBC       = 'Dirichlet';
expectedBC      = {'Dirichlet','Neumann'};
defaultGhostSize    = 0;

% check the order
if order < 1
    error('Error: polynomial order should be at least 1.')
end

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
addRequired(p,'order',validScalarPosNum);
% validate domain
validDomain = @(x) length(x)==2 && x(2)>x(1);
addParameter(p,'domain',defaultDomain,validDomain);
% validate type
addParameter(p,'type',defaultType,@(x) any(validatestring(x,expectedType)));
% number of lagrange elements
addParameter(p,'lag_elems',defaultLagElems,@(x)isnumeric(x) && isscalar(x) && (x >= 1));
% alternatively take a grid
addParameter(p,'lag_grid',defaultLagGrid,@(x)isnumeric(x));
% lagrange nodes type
addParameter(p,'lag_nodes',defaultLagNodes,@(x)any(validatestring(x,expectedLagNodes)));
% lagrange boundary condition
addParameter(p,'bc',defaultBC,@(x)any(validatestring(x,expectedBC)));
% lagrange boundary ghost cell size
addParameter(p,'ghost_size',defaultGhostSize,validScalarPosNum);
%
p.KeepUnmatched = true;
parse(p,order,varargin{:});

options = p.Results;

if ismember('domain',cellstr(p.UsingDefaults))
    switch options.type
        case{'Laguerre','Hermite'}
            options.domain = [-1, 1];
    end
end

% check inputs
switch options.type
    case{'Lagrange'}
        % check the ghost cell
        %
        if isempty(options.lag_grid)
            % old:
            %lag_grid = linspace(options.domain(1)+options.ghost_size, options.domain(2)-options.ghost_size, options.lag_elems+1);
            %
            % new:
            lag_grid = linspace(options.domain(1), options.domain(2), options.lag_elems+1);
            if (lag_grid(1) + options.ghost_size) > lag_grid(2)
                error('Ghost cell too large')
            end
            lag_grid(1) = lag_grid(1) + options.ghost_size;
            lag_grid(end) = lag_grid(end) - options.ghost_size;
        else
            lag_grid = options.lag_grid;
            disp('should check if the input grid is consistent, not implemented')
        end
        if options.ghost_size < eps
            options.bc  = 'none';
        end
        local       = setup_ref_lagrange(options.lag_nodes, order+1);
        def         = setup_lagrange(lag_grid, local, options.bc, options.ghost_size);
        def.local   = local;
        def.gs      = options.ghost_size;
        def.bc      = options.bc;
        def.order   = order;
        def.type    = options.type;
        
    otherwise
        def = setup_orthogonal(options.type, order);
        [def.nodes,J]   = reference2domain(def.ref_nodes, options.domain);
        def.omegas      = def.omegas/J;
end
def.domain = options.domain;

end