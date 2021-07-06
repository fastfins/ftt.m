classdef onedCDF
    
    properties
        tol
    end
    
    methods (Abstract)
        invert_cdf(obj)
        eval_cdf(obj)
        eval_cdf_deri(obj)
    end
        
    methods
        function obj = onedCDF(varargin)
            defaultErrTol = 1E-10;
            p = inputParser;
            %
            addOptional(p,'err_tol',defaultErrTol, @(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1));
            parse(p,varargin{:});
            obj.tol = p.Results.err_tol;
        end
    end
end