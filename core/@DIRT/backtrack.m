function gz = backtrack(Juz, Jux, gx)
% Given the deep Rosenblatt transport X = T(Z), where Z is the reference 
% random variable and X is the target random variable. and the gradient of 
% f(X), evaluate the gradient of f(T(Z)). 
% 
%   gz = BACKTRACK(Juz, Jux, gx)
%
%   Juz - cell array of the Jacobian of the diagonal map, U = CDF(Z), at 
%         each layer, dimension of each: d x n
%   Jux - cell array of the Jacobian of the SRT (CDF map), U = R(X), at 
%         each layer, dimension of each: d x (dn)
%   gx  - gradient of f w.r.t. x
%   gz  - gradient of f w.r.t. z

gz = gx;
[d,n] = size(gz);
k = length(Juz);
for l = 1:k
    for i = 1:n
        ind = (i-1)*d + (1:d);
        gz(:,i) = Juz{l}(:,i) .* ( tril(Jux{l}(:,ind))' \ gz(:,i) );
    end
end

end