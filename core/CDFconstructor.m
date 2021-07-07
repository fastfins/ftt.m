function cdf = CDFconstructor(poly, varargin)
% Select the one dimensional CDF function for a given collocation basis. 
%   cdf = CDFconstructor(poly)

if isa(poly, 'Lagrange1')
    cdf = Lagrange1CDF(poly, varargin{:});
elseif isa(poly, 'Lagrangep')
    cdf = LagrangepCDF(poly, varargin{:});
elseif isa(poly, 'Chebyshev1st') || isa(poly, 'Chebyshev2nd') || ...
        isa(poly, 'Legendre') || isa(poly, 'Jacobi11')
    cdf = ChebyshevCDF(poly, varargin{:});
elseif isa(poly, 'Fourier')
    cdf = FourierCDF(poly, varargin{:});
elseif isa(poly, 'Laguerre') || isa(poly, 'Hermite')
    disp('not yet implemented')
end

end