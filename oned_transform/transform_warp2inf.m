function x = transform_warp2inf(u, a, s, type)
%transform a d dimensional bounded variable in [-a,a]^d to [-inf,inf]^d
%
%Inputs:
%  u: dxn in [-a,a]^d
%
%  a: dx1 vector or a scalar, this defines the domain of x
%
%  s: dx1 vector or a scalar, this defines the scale of the transformation
%     For Gaussian, this defines the standard deviation. By default, s = 1, 
%     with increasing s, the rate of decay of the tails decreases
%
%  type: the reference measure of the output variable, 
%     with increasing tail masses, we have 'normal', 'logistic', 'cauchy'
%     default is 'normal'
%
%Outputs:
%  x: dxn, variable in [-inf,inf]^d
%

z = transform_warp2uni(u, a, 'normal');
x = transform_uni2inf(z, s, type);

end