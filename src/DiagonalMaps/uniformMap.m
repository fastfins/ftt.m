classdef uniformMap < diagonalMap
    
    methods
        function [u,Juz] = eval_cdf(obj, z)
            Juz = ones(size(z)) ./ (max(obj.domain) - min(obj.domain));
            u = (z - min(obj.domain)) .* Juz;
        end
        
        function  z = eval_icdf(obj, u)
            z = u .* (max(obj.domain) - min(obj.domain)) + min(obj.domain);
        end
        
        function [lf,glf] = log_pdf(obj, z)
            lf = -log(max(obj.domain) - min(obj.domain)) .* ones(1, size(z,2)); 
            lf(any((z<min(obj.domain))|(z>max(obj.domain)),1)) = -inf;
            glf = zeros(1, size(z,2)); 
        end
        
        function obj = uniformMap(domain)
            if (nargin<1)||(isempty(domain))
                domain = [0,1];
            end
            obj@diagonalMap(domain);
        end
        
    end
    
end