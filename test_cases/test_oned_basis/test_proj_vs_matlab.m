function test_proj_vs_matlab(p, f, a, b, plot_a, plot_b)

xs = linspace(plot_a, plot_b, 1000)';
[A,~]  = eval_basis(p, xs);

figure
plot(xs, f(xs), 'linewidth', 2)
hold on
title(func2str(f))

coef = p.node2basis*f(p.nodes(:));
plot(xs, A(:, 1:p.order+1)*coef)

c = zeros(size(coef));
for i = 1:length(c)
    c(i) = integral(@(xs) func(f, p, i, xs), a, b);
end

plot(xs, A(:, 1:p.order+1)*c)
plot(p.nodes, f(p.nodes), 'o', 'markersize', 10)
legend('true', 'mine', 'matlab', 'quad')

disp([coef c])

end

function f = func(fp, poly, i, xs)

[A,w]  = eval_basis(poly, xs(:));
f = A(:,i).*w.*fp(xs(:));

f = f';

end