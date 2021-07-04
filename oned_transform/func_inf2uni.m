function f = func_inf2uni(func, u, s, type)
%transform the coordinate of a function from [-inf,inf]^d to [0,1]^d, the 
%  evaluation contains the change of variable
%
%Inputs:
%  func: an input function defined in [-inf,inf]^d
%
%  u: dxn, variable in [0,1]^d with a uniform reference
%
%  s: dx1 vector or a scalar, this defines the scale of the transformation
%     For Gaussian, this defines the standard deviation. By default, s = 1, 
%     with increasing s, the rate of decay of the tails decreases
%
%Outputs:
%  f: []xn
%

[x, dzdx] = transform_uni2inf(u, s, type);
f = func(x)./prod(dzdx,1);

end