classdef uniformMap < diagonalMap
    
    methods
        function [u,Juz] = eval_cdf(obj, z)
            u = z;
            Juz = ones(size(z));
        end
        
        function  z = eval_icdf(obj, u)
            z = u;
        end
        
        function [lf,glf] = log_pdf(obj, z)
            lf = zeros(1, size(z,2)); 
            glf = zeros(1, size(z,2)); 
        end
        
        function obj = uniformMap()
            obj@diagonalMap([0,1]);
        end
    end
    
end