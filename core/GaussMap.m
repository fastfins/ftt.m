classdef GaussMap < diagonalMap
    
    methods
        function [u,Juz] = eval_cdf(obj, z)
            u = 0.5*( 1 + erf(z/sqrt(2)) ); 
            Juz = exp( -0.5*z.^2-0.5*log(2*pi) );
        end
        
        function  z = eval_icdf(obj, u)
            z = erfinv(u*2 - 1)*sqrt(2);
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
        end
        
        function [lf, glf] = log_pdf(obj, z)
            lf  = - 0.5*sum(z.^2, 1) - 0.5*log(2*pi)*size(z,1);
            glf = - z; 
        end
        
        function obj = GaussMap(domain)
            obj@diagonalMap(domain);
        end
    end
    
end