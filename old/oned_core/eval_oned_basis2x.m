function [f, w] = eval_oned_basis2x(def, coeff, x)
%Evaluate one dimensional basis functions (m dim.) on inputs x
%
%This function is not currently in use
%
%Inputs:
%  def:     definition of the function approximation
%  coeff:   coefficients of the polynomial
%  x:       vector of input variables
%
%Output:
%   f:      vector of outputs variables, dimension is aligned with x
%   w:      vector of the weighting function values, dimension is aligned with x
%
%Tiangang Cui, August, 2019

switch def.type
    case{'Lagrange'}
        f = eval_lagrange(def, coeff, x);
        f = reshape(f,size(x));
        w = ones(size(x));
        
    otherwise
        [b,w] = eval_orthogonal_basis(def.type, def.domain, def.order, x);
        f = reshape(b*coeff,size(x));
        w = reshape(w,size(x));
end

end


