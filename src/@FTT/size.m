function [d,rs,ns] = size(obj)
% Find the ranks of the TT cores and the degrees of freedom of
% the approximation basis for each coordinate. 
%   [d,rs,ns] = SIZE(tt)
%
%   d   - dimension of the function
%   rs  - ranks, r(0) = 1 is not included
%   ns  - number of collocation points for each dimension

d   = obj.d;
rs  = rank(obj);
ns  = ones(1,d);
for k = 1:d
    ns(k) = obj.oneds{k}.num_nodes;
end
end