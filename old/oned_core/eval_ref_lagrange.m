function f = eval_ref_lagrange(nodes, omega, f_at_x, x)
%Evaluate the reference Lagrange polynomail at input x for given function 
%value at nodes. node and omega are Lagrange nods and barycentric weights
%
%Tiangang Cui, August, 2019

tau = eps;
m   = length(x);
n   = length(nodes);
f   = zeros(m,1);

outside = nodes(1)-tau >= x | nodes(n)+tau <= x;
inside  = ~outside;
if sum(outside)
    disp('warning: points outside of the domain')
    f(outside) = 0;
end

tmp_x   = x(inside);
diff    = repmat(tmp_x(:), 1, n) - repmat(nodes(:)', length(tmp_x), 1);

% stablise
diff(abs(diff)<tau) = tau;
tmp_m   = repmat(omega(:)', length(tmp_x), 1) ./ diff;

% evaluation of the internal interpolation
f(inside)   = sum(repmat(f_at_x(:)',length(tmp_x),1).*tmp_m, 2)./sum(tmp_m, 2);

f = reshape(f, size(x));

end