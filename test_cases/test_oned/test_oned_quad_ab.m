% test the quadrature
a = 1
b = 5
max_order = 20;
polys = cell(4, max_order);
for order = 1:max_order
    polys{1, order} = Chebyshev1st(order, [a, b]);
    polys{2, order} = Chebyshev2nd(order, [a, b]);
    polys{3, order} = Legendre(order, [a, b]);
    polys{4, order} = Jacobi11(order, [a, b]);
end



%%%%%%%% case 1
test_bounded_ab(polys, @(y) sin(y*6*pi), a, b)

test_bounded_ab(polys, @(y) sin(y*5*pi), a, b)

test_bounded_ab(polys, @(y) sin(y*4.5*pi), a, b)

test_bounded_ab(polys, @(y) sqrt(y), a, b)

test_bounded_ab(polys, @(y) 1./sqrt(y), a, b)

test_bounded_ab(polys, @(y) 1./(2+y.^2), a, b);

