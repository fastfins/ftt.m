classdef ChebyshevCDF < Chebyshev2nd & SpectralCDF
        
    methods
        function obj = ChebyshevCDF(poly, varargin)
            obj@Chebyshev2nd(poly.order*2, poly.domain); 
            obj@SpectralCDF(varargin{:});
        end
    end
    
end