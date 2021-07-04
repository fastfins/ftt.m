function [def, options] = setup_oned_cdf(oned, varargin)
%
%Define data structure for building cdf functions from either the approximation
%of pdf (no guarantee of monotonicity) or the squared root of the pdf (guarantee
%of monotonicity).
%
%We transform the input polynomial for defining pdf to a Chebyshev polynomial
%(default is the 2nd), then the cdf can be obtained by indefinite integral
%of Chebyshev polynomial, which is analytically available.
%
%%%%%
%Input:
%
%oned:
%  The input oned polynomial
%
%%%%%
%Optional input parameters:
%
%squared:
%  indicate if the squared mode is used, default is true.
%
%cheby:
%  Type of Chebyshev polynomial used for building the cdf, choices are
%  'Chebyshev1st' and 'Chebyshev2nd', the default is 'Chebyshev2nd'
%
%err_tol:
%  Tolerance for solving the inverse transform, default is 1E-10
%
%method:
%  Root finding method for solving the inverse transform, choices are
%  'Regula_Falsi','Illinois'. The defauls is 'Regula_Falsi'
%
%%%%%
%Output:
%
%def:
%  A data structure contains:
%
%  nodes:       nodes for evaluating the pdf function
%  cdf_basis2node:
%               cdf basis functions defined at nodes for fast construction
%               of cdf grid.
%  squared:     indicate if the squared mode is used, default is true.
%  cheby:       type of Chebyshev polynomial used for building the cdf
%  err_tol:     tolerance for solving the inverse transform
%  method:      root finding method for solving the inverse transform
%
%  other components of this data structure are not directly used
%
%Tiangang Cui, August, 2019


defaultSquared  = true;
defaultCheby    = 'Chebyshev2nd';
expectedCheby   = {'Chebyshev1st','Chebyshev2nd'};
defaultErrTol   = 1E-10;
defaultMethod   = 'Regula_Falsi';
expectedMethod  = {'Regula_Falsi','Illinois'};

p = inputParser;
%
addRequired(p,'oned');
addParameter(p,'squared',defaultSquared,@(x) islogical(x) && isscalar(x));
addParameter(p,'cheby',  defaultCheby,  @(x) any(validatestring(x,expectedCheby)));
addParameter(p,'err_tol',defaultErrTol, @(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1));
addParameter(p,'method', defaultMethod, @(x) any(validatestring(x,expectedMethod)));

parse(p,oned,varargin{:});
options = p.Results;

if options.squared
    switch oned.type
        case{'Lagrange'}
            n_nodes = oned.order*2+1;
            local   = lag2cheby(options.cheby, oned.local.domain, n_nodes);
            def     = setup_lagrange_cdf(oned, local);
            def.local = local;
            def.order = oned.order*2;
            def.type  = oned.type;
            def.cheby = options.cheby;
            
        case{'Chebyshev1st', 'Chebyshev2nd', 'Legendre', 'Jacobi11'}
            order = oned.order*2;
            def   = setup_oned(order, 'type', options.cheby, 'domain', oned.domain);
            
        case{'Fourier'}
            order = oned.order*2+1;
            def   = setup_oned(order, 'type', 'Fourier', 'domain', oned.domain);
            
        case{'Laguerre', 'Hermite'}
            disp('not yet implemented')
    end
else
    switch oned.type
        case{'Lagrange'}
            n_nodes = oned.order+1;
            local   = lag2cheby(options.cheby, oned.local.domain, n_nodes);
            def     = setup_lagrange_cdf(oned, local);
            def.local = local;
            def.order = oned.order;
            def.type  = oned.type;
            def.cheby = options.cheby;
            
        case{'Legendre', 'Jacobi11'}
            def = setup_oned(oned.order, 'type', options.cheby, 'domain', oned.domain);
            
        case{'Laguerre', 'Hermite'}
            disp('not yet implemented')
        otherwise
            % case{'Chebyshev1st', 'Chebyshev2nd'}
            %   disp('input is already a chebyshev polynomial, no change of order needed')
            def = oned;
            
    end
end

% evaluate the intgrated basis function at nodes
switch oned.type
    case{'Lagrange'}
        def.cdf_basis2node = eval_oned_int_basis(options.cheby, def.local.domain, def.order, def.local.nodes);
        def.grid  = [def.elem_left, def.elem_right(end)];
    case{'Chebyshev1st', 'Chebyshev2nd', 'Legendre', 'Jacobi11'}
        def.sampling_nodes      = linspace(def.domain(1), def.domain(2), max(def.num_nodes*2, 200));
        def.sampling_nodes(1)   = def.domain(1)-eps;
        def.sampling_nodes(end) = def.domain(2)+eps;
        def.cdf_basis2node      = eval_oned_int_basis(options.cheby, def.domain, def.order, def.sampling_nodes);
    case{'Fourier'}
        def.sampling_nodes      = linspace(def.domain(1), def.domain(2), max(def.num_nodes*2, 200));
        def.sampling_nodes(1)   = def.domain(1)-eps;
        def.sampling_nodes(end) = def.domain(2)+eps;
        def.cdf_basis2node      = eval_oned_int_basis('Fourier', def.domain, def.order, def.sampling_nodes);
end

% copy properties over
def.squared     = options.squared;
def.err_tol     = options.err_tol;
def.method      = options.method;

end