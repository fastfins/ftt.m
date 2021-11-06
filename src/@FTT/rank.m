function rs = rank(obj)
% Find the ranks of the TT cores and the degrees of freedom of
% the approximation basis for each coordinate. 
%   rs = SIZE(tt)
%
%   rs  - ranks, r(0) = 1 is not included

rs  = ones(1,obj.d);
for k = 1:obj.d-1
    rs(k) = size(obj.cores{k+1}, 1);
end

end