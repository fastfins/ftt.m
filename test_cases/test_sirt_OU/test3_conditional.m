
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);

data = setup_ou_process(d, a);

func = @(x) eval_ou_process(data, x);
func2 = @(x1, x2, dir) eval_cond_ou_process(data, x1, x2, dir);

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
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);
opts{2} = FTToption('tt_method', 'random', 'sqrt_flag', true, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 5);

%%%%

for i = 1:4
    for j = 1:2
        tic;
        irts{i,j} = SIRT(func, d, polys{i}, opts{j}, 'debug_x', debug_x, 'sample_x', sample_x);
        toc
    end
end

% test m = 1, m = d-1, m = 8
% sample
z = rand(d, 1E4);
for i = 1:4
    for j = 1:2
        m = 8;
        ind1 = 1:m;
        ind2 = (m+1):d;
        xc1 = debug_x(ind1,:);
        xc2 = debug_x(ind2,:);
        % conditional sample
        zc1 = z(ind1, :);
        zc2 = z(ind2, :);
        %
        if irts{i,j}.int_dir ~= 1, irts{i,j} = marginalise(irts{i,j}, 1); end
        tic;[rc2,f] = eval_cirt(irts{i,j}, xc1, zc2);toc
        %
        figure;
        subplot(2,2,1);plot(abs(func2(xc1, rc2, 1) - f)); title('actual function value vs fft')
        subplot(2,2,2);plot(data.C(ind2, ind2) - cov(rc2')); title('actual covariance vs sample covariance')
        %
        if irts{i,j}.int_dir ~= -1, irts{i,j} = marginalise(irts{i,j}, -1); end
        tic;[rc1,f] = eval_cirt(irts{i,j}, xc2, zc1);toc
        %
        subplot(2,2,3);plot(abs(func2(rc1, xc2, -1) - f)); title('actual function value vs fft')
        subplot(2,2,4);plot(data.C(ind1, ind1) - cov(rc1')); title('actual covariance vs sample covariance')
    end
end

for i = 1:4
    for j = 1:2
        m = 8;
        indx = 1:m;
        indy = (m+1):d;
        x = debug_x(indx,1);
        my = data.C(indy,indx)*(data.C(indx,indx)\x); 
        Cy = data.C(indy,indy) - data.C(indy,indx)*(data.C(indx,indx)\data.C(indx,indy));
        % conditional sample
        zy = z(indy, :);
        %
        if irts{i,j}.int_dir ~= 1, irts{i,j} = marginalise(irts{i,j}, 1); end
        tic;[y,f] = eval_cirt(irts{i,j}, x, zy);toc
        %
        figure;
        subplot(2,3,1);plot(abs(func2(repmat(x,1,size(zy,2)), y, 1) - f)); title('actual function value vs fft')
        subplot(2,3,2);plot(Cy - cov(y')); title('actual covariance vs sample covariance')
        subplot(2,3,3);plot(my - mean(y,2)); title('actual mean vs sample mean')
        
        
        m = 8;
        indy = 1:m;
        indx = (m+1):d;
        x = debug_x(indx,1);
        my = data.C(indy,indx)*(data.C(indx,indx)\x); 
        Cy = data.C(indy,indy) - data.C(indy,indx)*(data.C(indx,indx)\data.C(indx,indy));
        % conditional sample
        zy = z(indy, :);
        %
        if irts{i,j}.int_dir ~= -1, irts{i,j} = marginalise(irts{i,j}, -1); end
        tic;[y,f] = eval_cirt(irts{i,j}, x, zy);toc
        %
        subplot(2,3,4);plot(abs(func2(y, repmat(x,1,size(zy,2)), -1) - f)); title('actual function value vs fft')
        subplot(2,3,5);plot(Cy - cov(y')); title('actual covariance vs sample covariance')
        subplot(2,3,6);plot(my - mean(y,2)); title('actual mean vs sample mean')
    end
end
