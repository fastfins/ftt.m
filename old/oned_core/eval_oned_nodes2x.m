function f = eval_oned_nodes2x(oned, x, tmp)

% x : row vector, 1 x n
% tmp: function at nodes, n_nodes x ??

switch oned.type
    case{'Chebyshev1st','Chebyshev2nd','Legendre','Jacobi11','Laguerre','Hermite','Fourier'}
        tmp = oned.node2basis*tmp;
end

f   = zeros(length(x), size(tmp,2));
mid = (x >= oned.domain(1)) & (x <= oned.domain(2));

if sum(mid) > 0
    b = eval_oned_basis(oned, x(mid));
    f(mid,:) = b*tmp;
end

end
