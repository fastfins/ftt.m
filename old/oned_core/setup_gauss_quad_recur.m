function [q_pts,q_weights] = setup_gauss_quad_recur(type, order)
% setup quadrature rules
%
% Inputs:
%   type:   type of orthogonal polynomial
%   order:  the order of polynomial, num. quad. points = order + 1
%
% Output:
%   q_pts:          quadrature points
%   quad_weights:   quadrature weights

% get the coefficients

k = (0:order)';

switch type
    
    case{'Chebyshev2nd'}
        recur.a = 2*ones(size(k));
        recur.b = zeros(size(k));
        recur.c = ones(size(k));
        
    case{'Legendre'}
        recur.a = (2*k+1)./(k+1);
        recur.b = zeros(size(k));
        recur.c = k./(k+1);
        
    case{'Jacobi11'}
        % the case alpha = beta = 1, for generating Lobatto points
        recur.a = (2*k+3).*(k+2)./(k+1)./(k+3);
        recur.b = zeros(size(k));
        recur.c = (k+2)./(k+3);
        
    case{'Laguerre'}
        recur.a = -1./(k+1);
        recur.b = (2*k+1)./(k+1);
        recur.c = k./(k+1);
        
    case{'Hermite'}
        recur.a = ones(size(k));
        recur.b = zeros(size(k));
        recur.c = k;
        
end

a = vpa(recur.a);
b = vpa(recur.b);
c = vpa(recur.c);

% assemble J matrix
% diagonal term
T0  = -b./a;
%{
Tp  = 1./a(1:end-1);
Tm  = c(2:end)./a(2:end);
T   = diag(T0,0) + diag(Tm,-1) + diag(Tp,1);
D   = diag(sqrt([1; a(2:end)./a(1)./cumprod(c(2:end))]));
%}
J1  = sqrt( c(2:end)./(a(1:end-1).*a(2:end)));
% J

%syms J
%J   = vpa( diag(T0,0) + diag(J1,1) + diag(J1,-1) );
J   = diag(T0,0) + diag(J1,1) + diag(J1,-1);

% off diagonal

[V,L]       = eig(vpa(J));
[q_pts,ind] = sort(diag(L));
q_pts       = double(q_pts(:));
q_weights   = double(V(1,ind)'.^2);

switch type
    case{'Chebyshev2nd'}
        q_weights   = q_weights*pi/2;
    case{'Jacobi11'}
        q_weights   = q_weights*4/3;
end

end

