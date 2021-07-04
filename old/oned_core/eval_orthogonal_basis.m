function [b, w] = eval_orthogonal_basis(type, domain, order, x)
% Evaluate one dimensional basis functions (m dim.) on inputs x
%
% Inputs:
%   def:    definition of the function approximation
%   x:      n x 1 vector of input variables
%
% Output:
%   b:      n x m, vector of outputs
%   w:      n x 1, vector of the weighting function values
%
%Tiangang Cui, August, 2019

switch type
    case{'Laguerre','Hermite'}
        J = 1;
    otherwise
        [x,J] = domain2reference(x,domain);
end

switch type
    case{'Chebyshev1st'}
        [b,w]   = eval_cheby(1, order, x);
    case{'Chebyshev2nd'}
        [b,w]   = eval_cheby(2, order, x);
        %{
        [b1,w1] = eval_oned_poly_recur(type, order, x);
        norm(b1(:)-b(:))
        norm(w1(:)-w(:))
        %}
    case{'Legendre','Jacobi11','Laguerre','Hermite'}
        [b,w]   = eval_oned_poly_recur(type, order, x);
        
    case {'Fourier'}
        m   = order+1;
        n   = m*2;
        is  = 2:m;
        ic  = (m+1):(n-1);
        
        b = zeros(length(x), n);
        b(:,1)  = ones(size(x(:)))*sqrt(0.5);
        b(:,is) = sin( x(:)*(1:order)*pi );
        b(:,ic) = cos( x(:)*(1:order)*pi );
        b(:,n)  = cos( x(:)*(m*pi) )*sqrt(0.5);
        
        w = ones(size(x));
        
end
w = w./J;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,w] = eval_cheby(type, order, x)
%
% Evaluate Chebyshev polynomials of the first and second kinds,
% for all input x, up to order n
%
% Inputs:
% order:
% x:    n_pts
%
% Output:
% f:    function outputs for each order at each x, n_pts x (n+1)
% w:    weight function at each x, n_pts x 1
%
% x = cos(theta), or theta = acos(x), x in [-1, 1]
% T_n(x) = cos( n * theta )
% w(x) = 1 / sqrt(1 - x^2)
%
% U_n(theta) = cos( (n+1) * theta )  / cos(theta)
% w(x) = sqrt(1 - x^2)
%

theta = real(acos(x(:)));
n     = reshape(0:order,1,[]);
switch type
    case{1}
        normalising = reshape([1,sqrt(2)*ones(1,order)]/sqrt(pi), 1,[]);
        f = cos(theta.*n) .* normalising;
    case{2}
        normalising = sqrt(2/pi);
        % deal with end points
        f = sin(theta.*(n+1)) ./ (sin(theta)/normalising);
        
        mask = abs(x+1) < eps;
        if sum(mask) > 0
            f(mask,:) = ((n+1).*(-1).^n)*normalising;
        end
        
        mask = abs(x-1) < eps;
        if sum(mask) > 0
            f(mask,:) = (n+1)*normalising;
        end
end

if nargout > 1
    switch type
        case{1}
            w = 1 ./ sqrt( 1 - x(:).^2 );
        case{2}
            w = sqrt( 1 - x(:).^2 );
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, w] = eval_oned_poly_recur(type, order, x)
% alpha = 1, beta = 1
% change of variable, x in [0, 1], then y = 2*x - 1
% all the normalising terms and recurrence terms can be pre-computed
k = 0:order;
switch type
    
    case{'Chebyshev2nd'}
        w = sqrt( 1 - x(:).^2 );
        normalising = sqrt(2/pi);
    case{'Legendre'}
        w = 0.5*ones(size(x));
        normalising = repmat( sqrt(2*(0:order)+1) , length(x), 1);
    case{'Jacobi11'}
        % the case alpha = beta = 1, for generating Lobatto points
        w = 1-x(:).^2;
        normalising = repmat( sqrt( (2*k+3).*(k+2)./(8*(k+1)) ), length(x), 1);
    case{'Laguerre'}
        w   = exp(-x);
        normalising = 1;
    case{'Hermite'}
        w   = exp(-0.5*x.^2)/sqrt(2*pi);
        normalising = repmat( sqrt(1./cumprod([1, 1:order])) , length(x), 1);
end

if order == 0
    f   = ones(length(x), 1).*normalising;
    return;
end

switch type
    
    case{'Chebyshev2nd'}
        a = 2*ones(size(k));
        b = zeros(size(k));
        c = ones(size(k));
    case{'Legendre'}
        a = (2*k+1)./(k+1);
        b = zeros(size(k));
        c = k./(k+1);
    case{'Jacobi11'}
        % the case alpha = beta = 1, for generating Lobatto points
        a = (2*k+3).*(k+2)./(k+1)./(k+3);
        b = zeros(size(k));
        c = (k+2)./(k+3);
    case{'Laguerre'}
        a = -1./(k+1);
        b = (2*k+1)./(k+1);
        c = k./(k+1);
        w   = exp(-x);
    case{'Hermite'}
        a = ones(size(k));
        b = zeros(size(k));
        c = k;
end

f = zeros(length(x), order+1);

f(:,1) = 1;
f(:,2) = (a(1)*x(:)+b(1)).*f(:,1);

for j = 2:order
    f(:,j+1) = (a(j)*x(:) + b(j)).*f(:,j) - c(j)*f(:,j-1);
end

f = f.*normalising;

end

