classdef FourierCDF < Fourier & SpectralCDF
    
    methods
        function obj = FourierCDF(poly, varargin)
            obj@Fourier(poly.order*2+1, poly.domain); 
            obj@SpectralCDF(varargin{:});
        end
    end
    
end