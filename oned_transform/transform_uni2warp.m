function [x, dudx] = transform_uni2warp(u, a, type)
%transform a d dimensional bounded variable in [0,1]^d to [-a,a]^d
%
%Inputs:
%  u: dxn in [0,1]^d
%
%  a: dx1 vector or a scalar, this defines the domain of x
%
%  type: the reference measure of the output variable, 
%     with increasing tail masses, we have 'normal', 'logistic', 'cauchy'
%     default is 'normal'
%
%Outputs:
%  x: dxn, variable in [-a,a]^d
%
%  dudx: dxn, inverse Jacobian of the transformation
%

switch type
    case {'normal'}
        tail = 1/2 + erf(-a/sqrt(2))/2;
        z = 1 - tail*2;
        x    = erfinv(2*u.*z-1)*sqrt(2);
        dudx = exp(-0.5*x.^2)./(z*sqrt(2*pi));
    case {'logistic'}
        tail = 1./(1+exp(a));
        z = 1 - tail*2;
        t = 1./(u.*z)-1;
        x = -log(t);
        dudx = (t./z)./(1+t).^2;
    case {'cauchy'}
        tail = 1/2 + atan(-a)/pi;
        z = 1 - tail*2;
        x = tan( (u.*z - 1/2)*pi );
        dudx = 1/( (1+x.^2).*(z*pi) );
end

end
