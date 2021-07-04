function [x, dudx] = transform_uni2inf(u, s, type)
%transform a d dimensional bounded variable in [0,1]^d to [-inf,inf]^d
%
%Inputs:
%  u: dxn in [0,1]^d
%
%  s: dx1 vector or a scalar, this defines the scale of the transformation
%     For normal, this defines the standard deviation. By default, s = 1, 
%     with increasing s, the rate of decay of the tails decreases
%
%  type: the reference measure of the output variable, 
%     with increasing tail masses, we have 'normal', 'logistic', 'cauchy'
%
%Outputs:
%  x: dxn, variable in [0,1]^d
%
%  dudx: dxn, inverse Jacobian of the transformation
%

switch type
    case {'normal'}
        t = erfinv(2*u-1);
        x = t*(sqrt(2)*s);
        dudx = exp(-t.^2)./(sqrt(2*pi)*s);
    case {'logistic'}
        t = 1./u-1;
        t(t<=0) = 1E-50;
        x = -log(t).*s;
        dudx = (t./s)./(1+t).^2;
    case {'cauchy'}
        t = tan( (u-1/2)*pi ); % t = x/s
        x = t.*s;
        dudx = 1./((1+t.^2).*(s*pi));
end

end