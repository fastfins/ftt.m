function fx = eval(obj, x, varargin)
% Evaluate the ftt function. The output is horizontally aligned.
%
%   x  - input variables, d x n.
%
%   pt - (optional) permuation index, where fTT is built for f(y),
%        y is the permutation of x s.t. x = y(p) and y = x(pt)
%        for a given p, pt can be obtained by pt(p) = 1:length(p);
%
%   fx - function values at x, mxn
%

% permute the input
if length(varargin) > 1
    error('too many optional arguments')
elseif length(varargin) == 1
    pt = varargin{1};
    x = x(pt,:);
end
fx  = eval_block(obj, x, obj.direction);
if obj.direction > 0
    fx  = fx';
end
if obj.opt.sqrt_flag
    fx  = fx.^2;
end
end