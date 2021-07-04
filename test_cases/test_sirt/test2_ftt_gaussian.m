
% setup the OU process
d = 20;
a = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = setup_ou_process(d, a);

func1 = @(x) eval_ou_process(data, x);

%%%%

debug_size = 1E4;
debug_x = data.B\randn(d, debug_size);
sample_x = data.B\randn(d, 1E3);

%%%%

% setup the reference polynomial
poly1 = setup_oned(21, 'type', 'Legendre', 'domain', [-5,5]);
poly2 = setup_oned(10, 'type', 'Fourier',  'domain', [-5,5]);
poly3 = setup_oned(5, 'type', 'Lagrange', 'lag_elems', 4, 'ghost_size', 1E-5, 'domain', [-5,5]);

opt1 = ftt_options('method', 'AMEN', 'ng_flag', true, 'oned_ref', poly2, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
%opt1 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly1, ...
%    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt1 = build_ftt(func1, d, [], opt1, 'debug_x', debug_x, 'sample_x', sample_x);

% squared version using AMEN
% AMEN has higher rank than random sampling with the same accuracy
opt2 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly2, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt2 = build_ftt(func1, d, [], opt2, 'debug_x', debug_x, 'sample_x', sample_x);
   

opt3 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly3, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt3 = build_ftt(func1, d, [], opt3, 'debug_x', debug_x, 'sample_x', sample_x);


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

