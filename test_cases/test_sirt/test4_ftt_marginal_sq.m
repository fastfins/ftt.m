
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);

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

%{
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
%}

ftt2 = ftt1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test 1
ind  = 1;
mind = 2:d;
fttm = build_ftt_marginal(ftt2, mind);
ftt_irt = build_irt_marginal(fttm);
tic;[r,f] = eval_irt_marginal(ftt_irt, rand(length(ind), 1E4));toc
fe   = eval_ou_process_marginal(data, ind, r);
figure; 
subplot(1,2,1);plot(abs(fe - f)/max(abs(fe)), '.'); 
subplot(1,2,2);plot(fe , f, '.'); 
title('actual function value vs fft sq')
disp(['actual covariance - sample covariance: ', num2str(data.C(ind, ind) - cov(r'))])

% test 2
ind  = d;
mind = 1:d-1;
fttm = build_ftt_marginal(ftt2, mind);
ftt_irt = build_irt_marginal(fttm);
tic;[r,f] = eval_irt_marginal(ftt_irt, rand(length(ind), 1E4));toc
fe   = eval_ou_process_marginal(data, ind, r);
figure; 
subplot(1,2,1);plot(abs(fe - f)/max(abs(fe)), '.'); 
subplot(1,2,2);plot(fe , f, '.'); 
title('actual function value vs fft sq')
disp(['actual covariance - sample covariance: ', num2str(data.C(ind, ind) - cov(r'))])

% test 3
ind  = 1:13;
mind = 14:d;
fttm = build_ftt_marginal(ftt2, mind);
ftt_irt = build_irt_marginal(fttm);
tic;[r,f] = eval_irt_marginal(ftt_irt, rand(length(ind), 1E4));toc
fe   = eval_ou_process_marginal(data, ind, r);
figure; 
subplot(1,2,1);plot(abs(fe - f)/max(abs(fe)), '.'); 
subplot(1,2,2);plot(fe , f, '.'); 
title('actual function value vs fft sq')
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')

% test 4
ind  = 13:d;
mind = 1:12;
fttm = build_ftt_marginal(ftt2, mind);
ftt_irt = build_irt_marginal(fttm);
tic;[r,f] = eval_irt_marginal(ftt_irt, rand(length(ind), 1E4));toc
fe   = eval_ou_process_marginal(data, ind, r);
figure; 
subplot(1,2,1);plot(abs(fe - f)/max(abs(fe)), '.'); 
subplot(1,2,2);plot(fe , f, '.'); 
title('actual function value vs fft sq')
figure; plot(data.C(ind, ind) - cov(r')); title('actual covariance vs sample covariance')

% test 5
ind  = [1:11, 16:d];
mind = 12:15;
fttm = build_ftt_marginal(ftt2, mind);
%ftt_irt = build_irt_marginal(fttm);
fe   = eval_ou_process_marginal(data, ind, debug_x(ind,:));
fa   = eval_ftt_marginal(fttm, debug_x(ind,:));
figure; 
subplot(1,2,1);plot(abs(fe - fa./data.norm)/max(abs(fe)), '.'); 
subplot(1,2,2);plot(fe , fa./data.norm, '.'); 
title('actual function value vs fft sq')

% test 6
ind  = [1:8, 16:d];
mind = 9:15;
fttm = build_ftt_marginal(ftt2, mind);
%ftt_irt = build_irt_marginal(fttm);
fe   = eval_ou_process_marginal(data, ind, debug_x(ind,:));
fa   = eval_ftt_marginal(fttm, debug_x(ind,:));
figure; 
subplot(1,2,1);plot(abs(fe - fa./data.norm)/max(abs(fe)), '.'); 
subplot(1,2,2);plot(fe , fa./data.norm, '.'); 
title('actual function value vs fft sq')

