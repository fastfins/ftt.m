function lag = setup_ref_lagrange(type, n)
%
%Define the reference Lagrange basis, in the reference domain [0,1]
%This function should not be explicitly used
% 
%%%%%
%Input: 
%
%type:     
%  Type of the interpolation points, options are 'Jacobi' and 'Chebyshev' (2nd)
%
%n:        
%  Number of interpolation points, should be greater than or equal to 2
%
%%%%%
%Output:   
%  lag, a data structure contains:
%
%  domain:     the reference domain
%  nodes:      all the interpolation points, 1xn
%  num_nodes:  number of nodes
%  omega:      barycentric weights of the Lagrange polynomial
%  mass:       reference mass matrix
%  weights:    reference weighting factors (integration of each basis 
%              function), 1xn 
%
%Tiangang Cui, August, 2019

if n < 2
    disp('We need more than two points to define Lagrange interpolation')
    return
end

lag.domain      = [0,1];
lag.nodes = zeros(1,n);
lag.nodes(1)    = 0;
lag.nodes(end)  = 1;
lag.num_nodes   = n;

if n > 2
    switch type
        case{'Jacobi'}
            def = setup_orthogonal('Jacobi11', n-3);
        case{'Chebyshev'}
            def = setup_orthogonal('Chebyshev2nd', n-3);
    end
    % in the interval [0, 1]
    lag.nodes(2:n-1) = 0.5*(def.ref_nodes+1);
end

% compute the local omega coefficients
lag.omega = zeros(1,n);
for j = 1:n
    ind     = true(n,1);
    ind(j)  = false;
    lag.omega(j) = 1./prod( lag.nodes(j) - lag.nodes(ind) );
end

% define the mass matrix
lag.mass     = build_local_mass(lag.nodes, lag.omega, lag.domain);

% setup the intergration of each basis
lag.weights  = build_local_weights(lag.nodes, lag.omega, lag.domain);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mass = build_local_mass(nodes, omega, domain)

n = length(nodes);
I = eye(n);

mass = zeros(n);
for i = 1:n
    for j = 1:n
        fij = @(x) eval_ref_lagrange(nodes, omega, I(:,i), x).*eval_ref_lagrange(nodes, omega, I(:,j), x);
        mass(i,j) = integral(fij, domain(1), domain(2));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weights = build_local_weights(nodes, omega, domain)

n = length(nodes);
I = eye(n);

weights = zeros(1,n);
for i = 1:n
    fi = @(x) eval_ref_lagrange(nodes, omega, I(:,i), x);
    weights(i) = integral(fi, domain(1), domain(2));
end

%{
% use build in quadrature, may suffer from accuracy
[q_nodes,q_weights,~] = setup_orthogonal('Legendre', 15);

weights2 = zeros(1,n);
for i = 1:n
    weights2(i) = sum(eval_ref_lagrange(nodes, omega, I(:,i), q_nodes).*q_weights);
end

norm(weights(:) - weights2(:))
%}

end