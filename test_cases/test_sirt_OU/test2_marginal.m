
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);

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

opts{1} = FTToption('tt_method', 'amen', 'max_als', 5, 'sqrt_flag', true, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);
opts{2} = FTToption('tt_method', 'random', 'max_als', 5, 'sqrt_flag', true, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);


for i = 1:4
    for j = 1:2
        tic;
        irts{i,j} = SIRT(func, d, polys{i}, opts{j}, 'debug_x', debug_x, 'sample_x', sample_x);
        toc
        tmp = round(irts{i,1}, 1E-2);
        irts{i,3} = SIRT(func, d, tmp, 'debug_x', debug_x);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should test ind = 1, ind = 1:(d-1) for > 0
% should test ind = d, ind = 2:d for < 0
% sample
z = rand(d, 1E4);
for i = 1:4
    for j = 1:3
        figure;
        % test 1
        ind  = 1:8;
        if irts{i,j}.int_dir ~= 1, irts{i,j} = marginalise(irts{i,j}, 1); end
        tic;[r,f] = eval_irt(irts{i,j}, z(ind,:));toc
        fx = eval_pdf(irts{i,j}, r);
        tic;z0 = eval_rt(irts{i,j}, r);toc
        norm(z(ind,:) - z0)
        norm(f - fx)
        fe   = eval_ou_process_marginal(data, ind, r);
        %
        subplot(2,3,1);plot(abs(fe - f)/max(abs(fe)), '.');
        subplot(2,3,2);plot(fe , f, '.');
        title('actual function value vs fft')
        subplot(2,3,3);plot(data.C(ind, ind) - cov(r'))
        
        % test 2
        ind  = d:-1:15;
        if irts{i,j}.int_dir ~= -1, irts{i,j} = marginalise(irts{i,j}, -1); end
        tic;[r,f] = eval_irt(irts{i,j}, z(ind,:));toc
        fx = eval_pdf(irts{i,j}, r);
        tic;z0 = eval_rt(irts{i,j}, r);toc
        norm(z(ind,:) - z0)
        norm(f - fx)
        %
        fe   = eval_ou_process_marginal(data, ind, r);
        subplot(2,3,4);plot(abs(fe - f)/max(abs(fe)), '.');
        subplot(2,3,5);plot(fe , f, '.');
        title('actual function value vs fft')
        subplot(2,3,6);plot(data.C(ind, ind) - cov(r'))
    end
end



