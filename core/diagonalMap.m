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
    
    methods
        u = eval_cdf(obj, z)
        % Map from reference to uniform
        
        z = eval_icdf(obj, u)
        % Map from uniform to reference
        
        logf = log_pdf(obj, z)
        % Return the log of the reference pdf
    end
    
end