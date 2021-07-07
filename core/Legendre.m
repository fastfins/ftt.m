classdef Legendre < recurr
    
    methods
        function obj = Legendre(order, varargin)
            [order,domain] = spectral.process_input(order,varargin{:});
            %
            k = double((0:order))';
            a = (2*k+1)./(k+1);
            b = zeros(size(k));
            c = k./(k+1);        
            normalising = reshape( sqrt(double(2*(0:order)+1)), 1, []);
            obj@recurr(order, domain, a, b, c, normalising);
        end
        
        function [f,w] = eval_ref_basis(obj, x)
            f = eval_ref_basis@recurr(obj, x);
            w = 0.5*ones(size(x));
        end
        
        function [f,w] = eval_ref_basis_deri(obj, x)
            f = eval_ref_basis_deri@recurr(obj, x);
            w = 0.5*ones(size(x));
        end
    end
end