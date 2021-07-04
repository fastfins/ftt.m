function [b, w] = eval_orthogonal_basis_deri(type, domain, order, x)
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
        [b,w]   = eval_cheby_deri(1, order, x);
    case{'Chebyshev2nd'}
        [b,w]   = eval_cheby_deri(2, order, x);
        %{
        [b1,w1] = eval_oned_poly_recur(type, order, x);
        norm(b1(:)-b(:))
        norm(w1(:)-w(:))
        %}
    case{'Legendre'}
        [b,w]   = eval_legendre_deri(order, x);
    case{'Jacobi11','Laguerre','Hermite'}
        disp('not implemented')
    case {'Fourier'}
        m   = order+1;
        n   = m*2;
        is  = 2:m;
        ic  = (m+1):(n-1);
        
        b = zeros(length(x), n);
        % orginal function
        %b(:,1)  = ones(size(x(:)))*sqrt(0.5);
        %b(:,is) = sin( x(:)*(1:order)*pi );
        %b(:,ic) = cos( x(:)*(1:order)*pi );
        %b(:,n)  = cos( x(:)*(m*pi) )*sqrt(0.5);
        
        %b(:,1)  = zeros(size(x(:)));
        c = reshape((1:order)*pi, 1, []);
        b(:,is) =  cos( x(:)*c ).*c;
        b(:,ic) = -sin( x(:)*c ).*c;
        b(:,n)  = -sin( x(:)*(m*pi) )*sqrt(0.5)*(m*pi);
        
        w = ones(size(x));
end
w = w./J;
b = b./J;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,w] = eval_cheby_deri(type, order, x)
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
% ( d theta / dx ) = inv(dx/dtheta) = inv( - sin(theta))
%
% sin(theta)^2 = (1-x^2)
%
% T_n(x) = cos( n * theta )
% w(x) = 1 / sqrt(1 - x^2)
% d T_n(x) / dx = n * sin( n * theta) / sin(theta)
%
% U_n(theta) = sin( (n+1) * theta )  / sin(theta)
% w(x) = sqrt(1 - x^2)
% d U_n(x) / dx = - (n+1) * cos((n+1)*theta) / sin(theta)^2 + sin((n+1)*theta) cos(theta) / sin(theta)^3
%               = ( - (n+1) * cos((n+1)*theta) + x sin((n+1)*theta) / sin(theta) ) / (1-x^2)
%               = ( (n+1) * cos((n+1)*theta) - x sin((n+1)*theta) / sin(theta) ) / (x^2-1)

theta = real(acos(x(:)));
n     = reshape(0:order,1,[]);
switch type
    case{1}
        normalising = reshape([1,sqrt(2)*ones(1,order)]/sqrt(pi), 1,[]);
        % original function: f = cos(theta.*n) .* normalising;
        
        % deal with end points
        mask = abs(x+1) < eps;
        if sum(mask) > 0
            theta(mask) = pi;
        end
        mask = abs(x-1) < eps;
        if sum(mask) > 0
            theta(mask) = 0;
        end
        f = (sin(theta*n).*n)./sin(theta);
        f = f.*normalising;
    case{2}
        normalising = sqrt(2/pi);
        % deal with end points
        % orginal function: f = sin(theta.*(n+1)) ./ (sin(theta)/normalising);
        
        % deal with end points
        mask = abs(x+1) < eps;
        if sum(mask) > 0
            theta(mask) = pi;
        end
        mask = abs(x-1) < eps;
        if sum(mask) > 0
            theta(mask) = 0;
        end
        
        f = ( cos(theta*(n+1)).*(n+1) - sin(theta*(n+1)).*(x(:)./sin(theta)) )./(x(:).^2-1);
        f = f.*normalising;
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

function [df, w] = eval_legendre_deri(order, x)
% alpha = 1, beta = 1
% change of variable, x in [0, 1], then y = 2*x - 1
% all the normalising terms and recurrence terms can be pre-computed
k = 0:order;
w = 0.5*ones(size(x));
normalising = repmat( sqrt(2*(0:order)+1) , length(x), 1);

if order == 0
    df  = zeros(length(x), 1).*normalising;
    return;
end

a = (2*k+1)./(k+1);
b = zeros(size(k));
c = k./(k+1);


f = zeros(length(x), order);
f(:,1) = 1;
f(:,2) = (a(1)*x(:)+b(1)).*f(:,1);

for j = 2:order
    f(:,j+1) = (a(j)*x(:) + b(j)).*f(:,j) - c(j)*f(:,j-1);
end

df = zeros(length(x), order+1);
%df(:,1) = 0;
%
for j = 1:order
    df(:,j+1) = j*f(:,j) + x(:).*df(:,j);
end

df = df.*normalising;

end

