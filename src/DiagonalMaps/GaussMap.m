classdef GaussMap < DiagonalMap
    
    properties
        left 
        right
    end
    
    methods
        function [u,Juz] = eval_cdf(obj, z)
            %%{
            u = 0.5*( 1 + erf(z/sqrt(2)) );
            %
            u(u > 1-eps) = 1-eps;
            u(u < eps) = eps;
            %
            Juz = exp( -0.5*z.^2-0.5*log(2*pi) );
            %%}
            %{
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
            u = 0.5*( 1 + erf(z/sqrt(2)) );
            u = (u-obj.left)/(obj.right-obj.left);
            %
            u(u > 1-eps) = 1-eps;
            u(u < eps) = eps;
            %
            Juz = exp( -0.5*z.^2-0.5*log(2*pi) )/(obj.right-obj.left);
            %}
        end
        
        function  z = eval_icdf(obj, u)
            %u = u*(obj.right-obj.left) + obj.left;
            z = erfinv(u*2 - 1)*sqrt(2);
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
        end
        
        function [lf, glf] = log_pdf(obj, z)
            s = size(z,1);
            %lf  = - 0.5*sum(z.^2, 1) - 0.5*log(2*pi)*s;
            lf  = - 0.5*sum(z.^2, 1) - 0.5*log(2*pi)*s - s*log((obj.right-obj.left));
            glf = - z; 
        end
        
        function z = random(obj, d, n)
            % pseudo random samples
            u = rand(d, n);
            %
            u = u*(obj.right-obj.left) + obj.left;
            z = erfinv(u*2 - 1)*sqrt(2);
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
        end
        
        function z = sobol(obj, d, n)
            % QMC samples using Sobol sequence
            S = sobolset(d);
            u = net(S, 2^ceil(log2(n)))';
            %
            u = u*(obj.right-obj.left) + obj.left;
            z = erfinv(u*2 - 1)*sqrt(2);
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
        end
        
        function obj = GaussMap(domain)
            obj@DiagonalMap(domain);
            % quantiles of the truncation
            obj.left  = normcdf(obj.domain(1));
            obj.right = normcdf(obj.domain(2));
        end
    end
    
end