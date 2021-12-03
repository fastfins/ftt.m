classdef DiagonalMap
    % diagonalMap class 
    %               - Superclass for all diagonal maps
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
        
        z = random(obj, d, n)
        % pseudo random samples
        
        z = sobol(obj, d, n)
        % QMC samples using Sobol sequence
        
        function obj = DiagonalMap(domain)
            obj.domain  = domain;
        end
        
        function y = get_domain(obj)
            y = obj.domain;
        end
    end
    
end