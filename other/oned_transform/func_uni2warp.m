function f = func_uni2warp(func, u, a)
%transform the coordinate of a function from [0,1]^d to [-a,a]^d, the 
%  evaluation contains the change of variable
%
%Inputs:
%  func: an input function defined in [-a,a]^d
%
%  u: dxn, variable in [-a,a]^d with a truncated normal reference
%
%  a: dx1 vector or a scalar, this defines the domain of x
%
%Outputs:
%  f: []xn
%

[z, dzdu] = transform_warp2uni(u, a, 'normal');
f = func(z).*prod(dzdu,1);

end