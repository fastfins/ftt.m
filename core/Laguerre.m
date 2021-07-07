classdef Laguerre < recurr
    
    methods
        function obj = Laguerre(order, varargin)
            [order,domain] = spectral.process_input(order,varargin{:});
            k = double((0:order))';
            a = -1./(k+1);
            b = (2*k+1)./(k+1);
            c = k./(k+1);
            normalising = ones(1, order+1);
            obj@recurr(order, domain, a, b ,c, normalising);
            %
        end
        
        function [f,w] = eval_ref_basis(obj, x)
            f = eval_ref_basis@recurr(obj, x);
            w = exp(-x);
        end
        
        function [f,w] = eval_ref_basis_deri(obj, x)
            f = eval_ref_basis_deri@recurr(obj, x);
            w = exp(-x);
        end
    end
    
end