function [b, w] = eval_oned_basis_deri(def, x)
%Evaluate one dimensional basis functions (m dim.) on inputs x
%
%Inputs:
%  def:     definition of the function approximation
%  x:       n x 1 vector of input variables
%
%Output:
%  b:       n x m, vector of outputs
%  w:       n x 1, vector of the weighting function values
%
%Tiangang Cui, August, 2019

switch def.type
    case{'Lagrange'}
        b = eval_lagrange_basis_deri(def, x(:));
        w = ones(size(x(:)));
    otherwise
        [b,w] = eval_orthogonal_basis_deri(def.type, def.domain, def.order, x(:));
end

end

% Note: need a fast way to


