function test_hermite(polys, f)

a = integral(f, -inf, inf);
max_quad = size(polys, 2);

err = zeros(1, max_quad);
for nquad = 1:max_quad
    err(1, nquad) = abs( sum(f( polys{nquad}.nodes)./ polys{nquad}.omegas.*polys{nquad}.weights) - a );
end

figure
semilogy(err', '.')
legend('Hermite')
title(func2str(f))

end