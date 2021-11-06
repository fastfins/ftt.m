function [u, dudx] = transform_inf2uni(x, s, type)
%transform a d dimensional unbounded variable in [-inf,inf]^d to [0,1]^d
%
%Inputs:
%  x: dxn, variable in [0,1]^d
%
%  s: dx1 vector or a scalar, this defines the scale of the transformation
%     For Gaussian, this defines the standard deviation. By default, s = 1, 
%     with increasing s, the rate of decay of the tails decreases
%
%  type: the reference measure of the input variable, 
%     with increasing tail masses, we have 'normal', 'logistic', 'cauchy'
%
%Outputs:
%  u: dxn in [0,1]^d
%
%  dudx: dxn, Jacobian of the transformation
%

switch type
    case {'normal'}
        t = x./(sqrt(2)*s);
        u = ( 1+erf(t) )/2;
        dudx = exp(-t.^2)./(sqrt(2*pi)*s);
    case {'logistic'}
        t = exp(-x./s);
        u = 1./(1+t);
        dudx = (t./s)./(1+t).^2;
    case {'cauchy'}
        t = x./s;
        u = 1/2 + atan(t)/pi;
        dudx = 1/( (1+t.^2).*(s*pi) );
end

end