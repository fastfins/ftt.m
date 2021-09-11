
% setup the OU process
d = 20;
a = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = setup_ou_process(d, a);
func = @(x) eval_ou_process(data, x);

%%%%

debug_size = 1E4;
debug_x = data.B\randn(d, debug_size);
sample_x = data.B\randn(d, 1E3);

%%%%

% setup the reference polynomial
polys{1} = Legendre(40, [-5,5]);
polys{2} = Fourier(20, [-5,5]);
polys{3} = Lagrangep(5, 8, [-5,5], 'ghost_size', 1E-5);
polys{4} = Lagrange1(40, [-5,5], 'ghost_size', 1E-5);

opts{1} = FTToption('tt_method', 'amen', 'sqrt_flag', true, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);
opts{2} = FTToption('tt_method', 'random', 'sqrt_flag', true, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);

for i = 1:4
    for j = 1:2
        tic;
        irts{i,j} = SIRT(func, d, polys{i}, opts{j}, 'debug_x', debug_x, 'sample_x', sample_x);
        toc
    end
    tmp = round(irts{i,1}, 1E-2);
    irts{i,3} = SIRT(func, d, tmp, 'debug_x', debug_x, 'sample_x', sample_x);
end

for i = 1:4
    for j = 1:3
        tic;
        irts{i,j} = marginalise(irts{i,j}, -1);
        toc
    end
end

% sample
z = rand(d, 1E4);
for i = 1:4
    for j = 1:3
        figure
        tic;[r,f] = eval_irt(irts{i,j}, z);toc
        tic;z0 = eval_rt(irts{i,j}, r);toc
        norm(z-z0, 'fro')
        %
        subplot(2,2,1); plot(abs(func(r)/data.norm - f), '.'); title('actual function vs fft')
        subplot(2,2,2); plot(func(r)/data.norm, f, '.');
        subplot(2,2,3); plot(data.C - cov(r')); title('actual covariance vs sample covariance')
    end
end

