classdef GaussMap < diagonalMap
    
    methods
        function u = eval_cdf(obj, z)
            u = 0.5*( 1 + erf(z/sqrt(2)) ); 
        end
        
        function  z = eval_icdf(obj, u)
            z = erfinv(u*2 - 1)*sqrt(2);
        end
        
        function logf = log_pdf(obj, z)
            logf = - 0.5*sum(z.^2, 1) - 0.5*log(2*pi)*size(z,1);
        end
    end
    
end