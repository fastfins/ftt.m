classdef uniformMap < diagonalMap
    
    methods
        function u = eval_cdf(obj, z)
            u = z;
        end
        
        function  z = eval_icdf(obj, u)
            z = u;
        end
        
        function logf = log_pdf(obj, z)
            logf = zeros(1, size(z,2)); 
        end
    end
    
end