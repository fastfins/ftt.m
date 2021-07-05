function [d,rs,ns] = size(obj)
% Find the ranks of the TT cores and the degrees of freedom of
% the approximation basis for each coordinate. 
%   [d,rs,ns] = SIZE(tt)
%
%   d   - dimension of the function
%   rs  - ranks, r(0) = 1 is not included
%   ns  - number of collocation points for each dimension

d   = length(obj.cores);
rs  = ones(1,d);
ns  = ones(1,d);
for k = 1:d-1
    rs(k) = size(obj.cores{k+1}, 1);
    ns(k) = obj.oneds{k}.num_nodes;
end
ns(d) = obj.oneds{d}.num_nodes;
end