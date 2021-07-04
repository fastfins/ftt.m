function def = setup_orthogonal(type, order)
%Define the reference orthogonal basis, for Chebyshev, Jacobi, and Legendre,
%it's defiend in the reference domain [-1,1]
%
%All the basis functions are normalised so their norm is 1
%Chebyshev polynomials and quadrature weights are defined using sine and
%cosine functions, other functions are defined using the recursive formula
%
%This function should not be called explicitly.
%
%%%%%
%Inputs:
%
%type:
%  Type of the polynomial, options are 'Chebyshev1st', 'Chebyshev2nd',
%  'Legendre', 'Jacobi11', 'Laguerre' and 'Hermite'
%
%order:
%  Max order of the polynomials
%
%%%%%
%Output:
%  def, a data structure contains:
%
%domain:        the reference domain
%ref_nodes:     quadrature points in the reference domain, nx1
%nodes:         quadrature points in the original domain, nx1, this is set
%               up in the setup_oned function
%order:         max order of the polynomials
%num_nodes:     number of nodes, = order + 1
%omegas:        weighting function at nodes, nx1
%weights:       quadrature weights, nx1
%node2basis:    a matrix transfer the function value at nodes to coefficients
%               of basis polynomials, nxn
%basis2node:    basis polynomials evaluated at quadrature nodes, nxn
%type:          type of the polynomial
%
%Tiangang Cui, August, 2019

switch type
    case{'Chebyshev1st'}
        % Recurrence does not hold for j = 1
        n = order + 1;
        def.ref_nodes   = reshape( sort( double( cos( vpa(pi)*(2*(1:n)-1)/(2*n) ) ), 'ascend'), [], 1);
        def.weights     = double( ones(size(def.ref_nodes))*vpa(pi)/n );
        
    case{'Chebyshev2nd'}
        n = order + 1;
        def.ref_nodes   = reshape( sort( double( cos( vpa(pi)*(1:n)/(n+1) ) ), 'ascend'), [], 1);
        def.weights     = double( sin( vpa(pi)*(1:n)/(n+1)).^2*vpa(pi)/(n+1) );
        
    case{'Legendre', 'Jacobi11', 'Laguerre', 'Hermite'}
        [def.ref_nodes,def.weights,~] = setup_gauss_quad_recur(type, order);
        
    case{'Fourier'}
        n = order*2 + 2;
        def.ref_nodes   = reshape( sort( (2/n)*(1:n) - 1, 'ascend'), [], 1);
        def.weights     = ones(size(def.ref_nodes))*(2/n);
end

[def.basis2node,def.omegas] = eval_orthogonal_basis(type, [-1,1], order, def.ref_nodes);
def.omegas      = reshape(def.omegas, size(def.ref_nodes));
def.weights     = reshape(def.weights, size(def.ref_nodes));
def.node2basis  = def.basis2node'*diag(def.weights);
def.num_nodes   = length(def.ref_nodes);
def.order       = order;
def.type        = type;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [q_pts,q_weights,recur] = setup_gauss_quad_recur(type, order)
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

