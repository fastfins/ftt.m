classdef Hermite < recurr

    methods
        function obj = Hermite(order, varargin)
            [order,domain] = spectral.process_input(order,varargin{:});
            k = double((0:order))';
            a = ones(size(k));
            b = zeros(size(k));
            c = k;
            normalising = reshape( sqrt(1./cumprod(double([1, 1:order]))), 1, []);
            obj = obj@recurr(order, domain, a, b, c, normalising);
        end
        
        function [f,w] = eval_ref_basis(obj, x)
            f = eval_ref_basis@recurr(obj, x);
            w = exp(-0.5*x.^2)/sqrt(2*pi);
        end
        
        function [f,w] = eval_ref_basis_deri(obj, x)
            f = eval_ref_basis_deri@recurr(obj, x);
            w = exp(-0.5*x.^2)/sqrt(2*pi);
        end
    end
    
end