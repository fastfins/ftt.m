function [u, dudx] = transform_warp2uni(x, a, type)
%transform a d dimensional bounded variable in [-a,a]^d to [0,1]^d
%
%Inputs:
%  x: dxn in [-a,a]^d
%
%  a: dx1 vector or a scalar, this defines the domain of x
%
%  type: the reference measure of the input variable, 
%     with increasing tail masses, we have 'normal', 'logistic', 'cauchy'
%     default is 'normal'
%
%Outputs:
%  u: dxn, variable in [0,1]^d
%
%  dudx: dxn, Jacobian of the transformation
%


switch type
    case {'normal'}
        tail = ( 1+erf(-a/sqrt(2)) )/2;
        z = 1 - tail*2;
        u = ( 1+erf(x/sqrt(2)) )./(z*2);
        dudx = exp(-0.5*x.^2)./(z*sqrt(2*pi));
    case {'logistic'}
        tail = 1./(1+exp(a));
        z = 1 - tail*2;
        t = exp(-x);
        u = 1./((1+t).*z);
        dudx = (t./z)./(1+t).^2;
    case {'cauchy'}
        tail = 1/2 + atan(-a)/pi;
        z = 1 - tail*2;
        u = ( 1/2 + atan(x)/pi )./z;
        dudx = 1/( (1+x.^2).*(z*pi) );
end

end