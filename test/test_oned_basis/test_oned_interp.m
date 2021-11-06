% test the interpolation

% first test the global one
nums = 2:11;
lags = cell(1, length(nums));

for j = 1:length(nums)
    lags{1, j} = Lagrangep(nums(j), 1);
end

%%%%%%%% case 1

test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) sin(y*5*pi))

test_interp_conv(lags(1,[5, 7, 8, 10]), @(y) sin(y*5*pi))


test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) sqrt(y));

test_interp_conv(lags(1,[5, 7, 8, 10]), @(y) sqrt(y));


test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) 1./sqrt(y));

test_interp_conv(lags(1,[5, 7, 8, 10]), @(y)  1./sqrt(y));


test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) log(y));

test_interp_conv(lags(1,[5, 7, 8, 10]), @(y) log(y));

%%%%%%%% case 2

xs  = linspace(0,1,1000)';
figure
hold on
lag = lags{1,end};
I   = eye(lag.num_nodes);
lf  = zeros(size(xs));
for j = 1:lag.num_nodes
    fi  = eval(lag, I(:,j), xs);
    lf  = lf + abs(fi);
    plot(xs, fi)
    plot(lag.nodes, I(:,j), 'o')
end
plot(xs, lf, 'linewidth', 2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
def = Lagrangep(5, 5, [-1.5, 1.5], 'bc', 'Neumann', 'ghost_size', 0.1);
plot(eval_basis(def, def.nodes))

figure
xs = linspace(def.domain(1), def.domain(2), 1000);
bs = eval_basis(def, xs);
for i = 1:size(bs,2), plot(xs, bs(:,i)); set(gca,'ylim', [-1,1]); title(num2str(i)); pause(0.2); end

%%%

def = Lagrangep(10, 3, [-1.5, 1.5], 'bc', 'Neumann', 'ghost_size', 0.1);
plot(eval_basis(def, def.nodes))

figure
xs = linspace(def.domain(1), def.domain(2), 1000);
bs = eval_basis(def, xs);
for i = 1:size(bs,2), plot(xs, bs(:,i)); set(gca,'ylim', [-1,1]); title(num2str(i)); pause(0.2); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


def = Lagrangep(10, 5, [-1.5, 1.5], 'bc', 'Dirichlet', 'ghost_size', 0.1);

xs = linspace(def.domain(1), def.domain(2), 1000);
% test interpolation property
f = @(y) sin(y*5*pi)+1;
fi = eval(def, f(def.nodes), xs);
plot(xs, f(xs), xs, fi, def.nodes,f(def.nodes(:)), 'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



def = Lagrangep(10, 5, [-5, -2], 'bc', 'Neumann', 'ghost_size', 0.1);
xs = linspace(def.domain(1), def.domain(2), 1000);
f = @(y) 1./log(-y) + sin(y*5*pi);
fi = eval(def, f(def.nodes), xs);
plot(xs, f(xs), xs, fi, def.nodes,f(def.nodes(:)), 'o')

%%%%
tic;
bs = eval_basis(def, xs);
for i = 1:10000
    f1 = bs*rand(def.num_nodes,1);
end
toc
tic;
bs = sparse(eval_basis(def, xs));
for i = 1:10000
    f1 = bs*rand(def.num_nodes,1);
end
toc
%
tic;
for i = 1:10000
    bs = eval_basis(def, xs);
    f1 = bs*rand(def.num_nodes,1);
end
toc
tic;
for i = 1:10000
    f = eval(def, rand(def.num_nodes,1), xs);
end
toc

%%%% test integration

err = zeros(5,6);

f = @(y) 1./log(abs(y)+0.01) + sin(y*5*pi);
a = integral(f, 2,4);
for i = 1:5
    for j = 1:6
        def = Lagrangep(2+j, 1+i, [2, 4], 'bc', 'Neumann');
        err(i,j) = abs(def.int_W*f(def.nodes(:)) - a);
    end
end
semilogy(err')
