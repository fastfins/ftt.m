classdef ChebyshevCDF < Chebyshev2nd & spectralCDF
        
    methods
        function obj = ChebyshevCDF(poly, varargin)
            obj@Chebyshev2nd(poly.order*2, poly.domain); 
            obj@spectralCDF(varargin{:});
        end
    end
    
end