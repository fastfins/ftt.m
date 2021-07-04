function test_laguerre(polys, f)

a = integral(f, 0, inf);
max_quad = size(polys, 2);

err = zeros(1, max_quad);
for nquad = 1:max_quad
    err(1, nquad) = abs( sum(f( polys{nquad}.nodes)./ polys{nquad}.omegas.*polys{nquad}.weights) - a );
end

figure
semilogy(err', '.')
legend('Laguerre')
title(func2str(f))

end