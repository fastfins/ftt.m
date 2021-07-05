
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);

data = setup_ou_process(d, a);

func1 = @(x) eval_ou_process(data, x);
func2 = @(x1, x2, dir) eval_cond_ou_process(data, x1, x2, dir);

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
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);
opts{2} = FTToption('tt_method', 'random', 'sqrt_flag', true, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);

for i = 1:4
    for j = 1:2
        tic;
        irts{i,j} = SIRT(func, d, polys{i}, opts{j}, 'debug_x', debug_x, 'sample_x', sample_x);
        irts{i,j} = marginalise(irts{i,j}, 1);
        toc
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
m = 1;
xc1 = debug_x(1:m,:);
xc2 = debug_x((d-m+1):d,:);
% conditional sample
z = rand(d-m, 1E4);

%%%%

ftt_int1 = build_irt(ftt1, 'dir', -1);
tic;[r, f] = eval_cond_irt(ftt_int1, xc2, z);toc
%
figure; plot(abs(func2(r, xc2, -1)/data.norm - f)); title('actual function value vs fft')
ind = 1:(d-m);
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')


% dir > 0, order of samples (xc, r), conditional (r | xc)
ftt_int1 = build_irt(ftt1, 'dir', 1);
tic;[r, f] = eval_cond_irt(ftt_int1, xc1, z);toc
figure; plot(abs(func2(xc1, r, 1)/data.norm - f)); title('actual function value vs fft')
ind = (m+1):d;
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')

%%%%


ftt_int2 = build_irt(ftt2, 'dir', -1);
tic;[r, f] = eval_cond_irt(ftt_int2, xc2, z);toc
%
figure; plot(abs(func2(r, xc2, -1)/data.norm - f)); title('actual function value vs fft')
ind = 1:(d-m);
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')


% dir > 0, order of samples (xc, r), conditional (r | xc)
ftt_int2 = build_irt(ftt2, 'dir', 1);
tic;[r, f] = eval_cond_irt(ftt_int2, xc1, z);toc
figure; plot(abs(func2(xc1, r, 1)/data.norm - f)); title('actual function value vs fft')
ind = (m+1):d;
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')


%%%%
ftt_int3 = build_irt(ftt3, 'dir', -1);
tic;[r, f] = eval_cond_irt(ftt_int3, xc2, z);toc
%
figure; plot(abs(func2(r, xc2, -1)/data.norm - f)); title('actual function value vs fft sq')
ind = 1:(d-m);
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')

ftt_int3 = build_irt(ftt3, 'dir', 1);
tic;[r, f] = eval_cond_irt(ftt_int3, xc1, z);toc
figure; plot(abs(func2(xc1, r, 1)/data.norm - f)); title('actual function value vs fft sq')
ind = (m+1):d;
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')

%%%%%%%%%%%%%%%
