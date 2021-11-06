classdef Jacobi11 < Recurr
    
    methods
        function obj = Jacobi11(order, varargin)
            [order,domain] = Spectral.process_input(order,varargin{:});
            %
            k = double((0:order))';
            a = (2*k+3).*(k+2)./(k+1)./(k+3);
            b = zeros(size(k));
            c = (k+2)./(k+3);
            normalising = reshape( sqrt( (2*k+3).*(k+2)./(8*(k+1)) ), 1, []);
            obj@Recurr(order, domain, a, b, c, normalising);
            %
            obj.weights = obj.weights*4/3;
            obj.node2basis = obj.node2basis*4/3;
        end
        
        function [f,w] = eval_ref_basis(obj, x)
            f = eval_ref_basis@Recurr(obj, x);
            w = 1-x.^2;
        end
        
        function [f,w] = eval_ref_basis_deri(obj, x)
            f = eval_ref_basis_deri@Recurr(obj, x);
            w = 1-x.^2;
        end
    end
end