function cdf = CDFconstructor(poly, varargin)

switch get_name(poly)
    case{'Lagrange1'}
        cdf = Lagrange1CDF(poly, varargin{:});
        
    case{'Lagrangep'}
        cdf = LagrangepCDF(poly, varargin{:});
        
    case{'Chebyshev1st', 'Chebyshev2nd', 'Legendre', 'Jacobi11'}
        cdf = ChebyshevCDF(poly, varargin{:});
        
    case{'Fourier'}
        cdf = FourierCDF(poly, varargin{:});
        
    case{'Laguerre', 'Hermite'}
        disp('not yet implemented')
end

end