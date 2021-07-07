function test_bounded(polys, f)

a = integral(f, 0, 1);
max_quad = size(polys, 2);

err = zeros(4, max_quad);
for nquad = 1:max_quad
    err(1, nquad) = abs( sum(f( polys{1,nquad}.nodes)./ polys{1,nquad}.omegas.*polys{1,nquad}.weights) - a );
    err(2, nquad) = abs( sum(f( polys{2,nquad}.nodes)./ polys{2,nquad}.omegas.*polys{2,nquad}.weights) - a );
    err(3, nquad) = abs( sum(f( polys{3,nquad}.nodes)./ polys{3,nquad}.omegas.*polys{3,nquad}.weights) - a );
    err(4, nquad) = abs( sum(f( polys{4,nquad}.nodes)./ polys{4,nquad}.omegas.*polys{4,nquad}.weights) - a );
end

names = cell(1, 4);
for i = 1:4
    names{i} = class(polys{i,1});
end

figure
semilogy(err', '-o')
title(func2str(f))
legend(names)

end