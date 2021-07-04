% test the quadrature


max_order = 20;
polys = cell(6, max_order);
for order = 1:max_order
    polys{1, order} = Chebyshev1st(order);
    polys{2, order} = Chebyshev2nd(order);
    polys{3, order} = Legendre(order);
    polys{4, order} = Jacobi11(order);
    polys{5, order} = Hermite(order, [-1, 1]);
    polys{6, order} = Laguerre(order, [-1, 1]);
end



%%%%%%%% case 1
test_bounded(polys, @(y) sin(y*6*pi))

test_bounded(polys, @(y) sin(y*5*pi))

test_bounded(polys, @(y) sin(y*4.5*pi))

test_bounded(polys, @(y) sqrt(y));

test_bounded(polys, @(y) 1./sqrt(y));

test_bounded(polys, @(y) 1./(2+y.^2));

test_bounded(polys, @(y) 10*y.^10 );

test_bounded(polys,  @(y) y.^10.*sqrt(y - y.^2 ));

test_bounded(polys,  @(y) y.^10./sqrt(y - y.^2 ));

test_bounded(polys, @(y) log(y));

test_bounded(polys, @(y) log(y+1));

%%%%%%% case 2
test_hermite(polys(5,:), @(y) exp(-6* y.^2 -y.^4));

test_hermite(polys(5,:), @(y) exp(-y.^2 - 0.5*y.^4) .* exp( - 0.5*y.^2) );

%%%%%%% case 3
test_laguerre(polys(6,:), @(y) exp(-y.^2 - 0.5*y.^4) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% discrete orthogonality 
j = 5;
for i = 1:6
    disp(i)
    A = polys{i,j}.basis2node; 
    disp(A'*diag(polys{i,j}.weights)*A)
end

j = 20;
for i = 1:6
    disp(i)
    A = polys{i,j}.basis2node; 
    figure
    plot(A'*diag(polys{i,j}.weights)*A)
    disp(mean(abs(diag(A'*diag(polys{i,j}.weights)*A))))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% projection test
j = 5;
for i = 1:3
    test_proj_vs_matlab(polys{i,j}, @(y) log(y) , 0, 1, 0, 1 );
end

test_proj_vs_matlab(polys{5,5}, @(x) 1./(1+x.^2) , -inf, inf, -6, 6 );

