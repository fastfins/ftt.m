
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);

data = setup_ou_process(d, a);

func1 = @(x) eval_ou_process(data, x);
func2 = @(x1, x2, dir) eval_cond_ou_process(data, x1, x2, dir);

%%%%

debug_size = 1E4;
debug_x = data.B\randn(d, debug_size);
sample_x = data.B\randn(d, 1E3);

%%%%

poly1 = setup_oned(21, 'type', 'Legendre', 'domain', [-5,5]);
poly2 = setup_oned(10, 'type', 'Fourier',  'domain', [-5,5]);
poly3 = setup_oned(5, 'type', 'Lagrange', 'lag_elems', 4, 'ghost_size', 1E-5, 'domain', [-5,5]);

opt1 = ftt_options('method', 'AMEN', 'ng_flag', true, 'oned_ref', poly1, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt1 = build_ftt(func1, d, [], opt1, 'debug_x', debug_x);

% squared version using AMEN
% AMEN has higher rank than random sampling with the same accuracy
opt2 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly2, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt2 = build_ftt(func1, d, [], opt2, 'debug_x', debug_x, 'sample_x', sample_x);
   

opt3 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly3, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt3 = build_ftt(func1, d, [], opt3, 'debug_x', debug_x, 'sample_x', sample_x);


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
