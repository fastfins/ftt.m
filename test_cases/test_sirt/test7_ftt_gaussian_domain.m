
% setup the OU process
d = 20;
a = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = setup_ou_process(d, a);
tr = 3;

func1 = @(x) eval_ou_process(data, x);
func2 = @(u) func_inf2warp(func1, u, tr, 1, 'cauchy');
func3 = @(u) func_inf2warp(func1, u, tr, 1, 'logistic');
func4 = @(u) func_inf2warp(func1, u, tr, 0.5, 'logistic');
%%%%

debug_size = 1E4;
debug_x = data.B\randn(d, debug_size);
sample_x = data.B\randn(d, 1E3);

debug_x = tr*1.95*debug_x/(max(debug_x(:)) - min(debug_x(:)));
sample_x = tr*1.95*sample_x/(max(sample_x(:)) - min(sample_x(:)));

%%%%

% setup the reference polynomial
poly1 = setup_oned(5, 'type', 'Lagrange', 'lag_elems', 4, 'ghost_size', 1E-5, 'domain', [-tr,tr]);
poly2 = setup_oned(10, 'type', 'Fourier',  'domain', [-tr,tr]);

opt1 = ftt_options('method', 'AMEN', 'ng_flag', true, 'oned_ref', poly1, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 15, 'max_als', 15);
opt2 = ftt_options('method', 'AMEN', 'ng_flag', true, 'oned_ref', poly2, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 15, 'max_als', 15);


ftt11 = build_ftt(func1, d, [], opt1, 'debug_x', debug_x, 'sample_x', sample_x);
ftt12 = build_ftt(func2, d, [], opt1, 'debug_x', debug_x, 'sample_x', sample_x);
ftt13 = build_ftt(func3, d, [], opt1, 'debug_x', debug_x, 'sample_x', sample_x);
ftt14 = build_ftt(func4, d, [], opt1, 'debug_x', debug_x, 'sample_x', sample_x);
   

ftt21 = build_ftt(func1, d, [], opt2, 'debug_x', debug_x, 'sample_x', sample_x);
ftt22 = build_ftt(func2, d, [], opt2, 'debug_x', debug_x, 'sample_x', sample_x);
ftt23 = build_ftt(func3, d, [], opt2, 'debug_x', debug_x, 'sample_x', sample_x);
ftt24 = build_ftt(func4, d, [], opt2, 'debug_x', debug_x, 'sample_x', sample_x);

ftt_int11 = build_irt(ftt11);
ftt_int12 = build_irt(ftt12);
ftt_int13 = build_irt(ftt13);
ftt_int14 = build_irt(ftt14);

ftt_int21 = build_irt(ftt21);
ftt_int22 = build_irt(ftt22);
ftt_int23 = build_irt(ftt23);
ftt_int24 = build_irt(ftt24);

%{
% sample
z = rand(d, 1E4);
% sample, squared version, 4 seconds for 1E4 samples
figure;
ftt_int1 = build_irt(ftt1, 'dir', -1);
tic;[r, f] = eval_irt(ftt_int1, z);toc
tic;z0 = eval_rt(ftt_int1, r);toc
norm(z-z0, 'fro')
subplot(2,3,1); plot(abs(func1(r)/data.norm - f), '.'); title('actual function vs fft')
subplot(2,3,2); plot(func1(r)/data.norm, f, '.');
subplot(2,3,3); plot(data.C - cov(r')); title('actual covariance vs sample covariance')

ftt_int1 = build_irt(ftt1, 'dir', 1);
tic;[r, f] = eval_irt(ftt_int1, z);toc
tic;z0 = eval_rt(ftt_int1, r);toc
norm(z-z0, 'fro')
subplot(2,3,4); plot(abs(func1(r)/data.norm - f), '.'); 
subplot(2,3,5); plot(func1(r)/data.norm, f, '.');
subplot(2,3,6); plot(data.C - cov(r')); 

% sample, squared version, 4 seconds for 1E4 samples
figure
ftt_int2 = build_irt(ftt2, 'dir', -1);
tic;[r, f] = eval_irt(ftt_int2, z);toc
tic;z0 = eval_rt(ftt_int2, r);toc
norm(z-z0, 'fro')
subplot(2,3,1); plot(abs(func1(r)/data.norm - f), '.'); title('actual function vs fft')
subplot(2,3,2); plot(func1(r)/data.norm, f, '.');
subplot(2,3,3); plot(data.C - cov(r')); title('actual covariance vs sample covariance')

ftt_int2 = build_irt(ftt2, 'dir', 1);
tic;[r, f] = eval_irt(ftt_int2, z);toc
tic;z0 = eval_rt(ftt_int2, r);toc
norm(z-z0, 'fro')
subplot(2,3,4); plot(abs(func1(r)/data.norm - f), '.'); 
subplot(2,3,5); plot(func1(r)/data.norm, f, '.');
subplot(2,3,6); plot(data.C - cov(r')); 

figure
ftt_int3 = build_irt(ftt3, 'dir', -1);
tic;[r, f] = eval_irt(ftt_int3, z);toc
tic;z0 = eval_rt(ftt_int3, r);toc
norm(z-z0, 'fro')
subplot(2,3,1); plot(abs(func1(r)/data.norm - f), '.'); title('actual function vs fft')
subplot(2,3,2); plot(func1(r)/data.norm, f, '.');
subplot(2,3,3); plot(data.C - cov(r')); title('actual covariance vs sample covariance')

ftt_int3 = build_irt(ftt3, 'dir', 1);
tic;[r, f] = eval_irt(ftt_int3, z);toc
tic;z0 = eval_rt(ftt_int3, r);toc
norm(z-z0, 'fro')
subplot(2,3,4); plot(abs(func1(r)/data.norm - f), '.'); 
subplot(2,3,5); plot(func1(r)/data.norm, f, '.');
subplot(2,3,6); plot(data.C - cov(r')); 
%%%%%%%%%%%%%%%
%}
