classdef GaussMap < diagonalMap
    
    properties
        left 
        right
    end
    
    methods
        function [u,Juz] = eval_cdf(obj, z)
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
            u = 0.5*( 1 + erf(z/sqrt(2)) ); 
            u = (u-obj.left)/(obj.right-obj.left);
            %
            Juz = exp( -0.5*z.^2-0.5*log(2*pi) )/(obj.right-obj.left);
        end
        
        function  z = eval_icdf(obj, u)
            u = u*(obj.right-obj.left) + obj.left;
            z = erfinv(u*2 - 1)*sqrt(2);
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
        end
        
        function [lf, glf] = log_pdf(obj, z)
            s = size(z,1);
            lf  = - 0.5*sum(z.^2, 1) - 0.5*log(2*pi)*s - s*log((obj.right-obj.left));
            glf = - z; 
        end
        
        function obj = GaussMap(domain)
            obj@diagonalMap(domain);
            % quantiles of the truncation
            obj.left  = normcdf(obj.domain(1));
            obj.right = normcdf(obj.domain(2));
        end
    end
    
end