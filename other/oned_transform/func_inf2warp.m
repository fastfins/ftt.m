function f = func_inf2warp(func, u, a, s, type)
%transform the coordinate of a function from [-inf,inf]^d to [-a,a]^d, the 
%  evaluation contains the change of variable
%
%Inputs:
%  func: an input function defined in [-inf,inf]^d
%
%  u: dxn, variable in [-a,a]^d with a truncated normal reference
%
%  a: dx1 vector or a scalar, this defines the domain of x
%
%  s: dx1 vector or a scalar, this defines the scale of the transformation
%     For Gaussian, this defines the standard deviation. By default, s = 1, 
%     with increasing s, the rate of decay of the tails decreases
%
%Outputs:
%  f: []xn
%

[z, dzdu] = transform_warp2uni(u, a, 'normal');
[x, dzdx] = transform_uni2inf(z, s, type);

f = func(x).*prod(dzdu./dzdx,1);

end