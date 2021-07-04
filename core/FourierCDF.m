classdef FourierCDF < Fourier & spectralCDF
    
    methods
        function obj = FourierCDF(poly, varargin)
            obj@Fourier(poly.order*2+1, poly.domain); 
            obj@spectralCDF(varargin{:});
        end
    end
    
end