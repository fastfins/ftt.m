classdef diagonalMap
    % diagonalMap class 
    %               - Superclass for all one dimensional basis for
    %                 evaluating CDF and inverse of the CDF
    %
    % diagonalMap Methods:
    %   eval_cdf    - Map from reference to uniform
    %   eval_icdf   - Map from uniform to reference
    %   log_pdf     - Return the log of the reference pdf
    %
    % See also uniformMap and GaussMap
    
    properties
        domain
    end
    
    methods
        u = eval_cdf(obj, z)
        % Map from reference to uniform
        
        z = eval_icdf(obj, u)
        % Map from uniform to reference
        
        logf = log_pdf(obj, z)
        % Return the log of the reference pdf
        
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
        
        function obj = diagonalMap(domain)
            obj.domain  = domain;
        end
        
        function y = get_domain(obj)
            y = obj.domain;
        end
    end
    
end