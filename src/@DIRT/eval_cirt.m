function [y,logf] = eval_cirt(obj, x, zy)
% Evaluate the inverse of the conditional deep Rosenblatt transport 
% Y | X = R^{-1}(Z, X), where X is given, (X, Y) jointly follow the target 
% distribution represented by SIRT, Z is uniform. 
%   [y,logf] = EVAL_CIRT(irt, x, z)
%
%   x    - given variables following the marginal of the target distribution
%   z    - reference random variables, m x n
%   y    - random variable drawn from Y | X
%   logf - log pdf of Y | X
%

d  = length(obj.d);
ny = size(zy,2);
dy = size(zy,1);
nx = size(x,2);
dx = size(x,1);

if dx == 0 || dy == 0
    error('dimension of x or the dimension of z should be nonzero')
end

if dy + dx ~= d
    error('dimension of x and the dimension of z mismatch the dimension of the joint')
end

zx = eval_rt (obj, x);
px = eval_pdf(obj, x);

if nx ~= ny
    if nx == 1
        zx = repmat(zx, 1, ny);
    else
        error('number of x and the number of z mismatch')
    end
end

[tmp, lf] = eval_irt(obj, [zx; zy]);
y = tmp((1:dy)+dx,:);
logf = lf - log(px);

end
