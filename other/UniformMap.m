classdef UniformMap < DiagonalMap
    
    methods
        function [u,Juz] = eval_cdf(obj, z)
            Juz = ones(size(z)) ./ (max(obj.domain) - min(obj.domain));
            u = (z - min(obj.domain)) .* Juz;
        end
        
        function  z = eval_icdf(obj, u)
            u(u > 1-eps) = 1-eps;
            u(u < eps) = eps;
            z = u .* (max(obj.domain) - min(obj.domain)) + min(obj.domain);
        end
        
        function [lf,glf] = log_pdf(obj, z)
            lf = -log(max(obj.domain) - min(obj.domain)) .* ones(1, size(z,2)); 
            %lf(any((z<min(obj.domain))|(z>max(obj.domain)),1)) = -inf;
            glf = zeros(1, size(z,2)); 
        end
        
        function z = random(obj, d, n)
            % pseudo random samples
            u = rand(d, n);
            z = eval_icdf(obj, u);
        end
        
        function z = sobol(obj, d, n)
            % QMC samples using Sobol sequence
            S = sobolset(d);
            u = net(S, 2^ceil(log2(n)));
            z = eval_icdf(obj, u');
        end
        
        function obj = UniformMap(domain)
            if (nargin<1)||(isempty(domain))
                domain = [0,1];
            end
            obj@DiagonalMap(domain);
        end
        
    end
    
end