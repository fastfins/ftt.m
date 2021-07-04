
% define three pdfs
%

gauss   = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
laplace = @(x) exp(-abs(x))/2;
cauchy  = @(x) 1./( pi*(1+x.^2) );

xs = linspace(-5, 5, 1000);

gs = gauss  (xs);
ls = laplace(xs);
cs = cauchy (xs);

figure(1)
plot(xs, gs, xs, ls, xs, cs)
legend('normal', 'laplace', 'cauchy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure(2);
set(f, 'position', [100, 100, 800, 800])
s = [0.5, 1, 2];
us = linspace(1E-16,1-1E-16,1000);

for i = 1:length(s)
    gs = func_inf2uni(gauss,   us, s(i), 'cauchy');
    ls = func_inf2uni(laplace, us, s(i), 'cauchy');
    cs = func_inf2uni(cauchy,  us, s(i), 'cauchy');
    
    subplot(length(s),3,(i-1)*3+1)
    plot(us, gs)
    if i == 1
        title('gauss, using cauchy')
    end
    ylabel(['s=' num2str(s(i))])
    subplot(length(s),3,(i-1)*3+2)
    plot(us, ls)
    if i == 1
        title('laplace, using cauchy')
    end
    subplot(length(s),3,(i-1)*3+3)
    plot(us, cs)
    if i == 1
        title('cauchy, using cauchy')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 5;
f = figure(3);
set(f, 'position', [100, 100, 800, 800])
ys = linspace(-5, 5, 1000);

for i = 1:length(s)
    gs = func_inf2warp(gauss,   ys, a, s(i), 'cauchy');
    ls = func_inf2warp(laplace, ys, a, s(i), 'cauchy');
    cs = func_inf2warp(cauchy,  ys, a, s(i), 'cauchy');
    
    subplot(length(s),3,(i-1)*3+1)
    plot(ys, gs)
    if i == 1
        title('gauss, using cauchy')
    end
    ylabel(['s=' num2str(s(i))])
    subplot(length(s),3,(i-1)*3+2)
    plot(ys, ls)
    if i == 1
        title('laplace, using cauchy')
    end
    subplot(length(s),3,(i-1)*3+3)
    plot(ys, cs)
    if i == 1
        title('cauchy, using cauchy')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure(4);
set(f, 'position', [100, 100, 800, 800])
s = [0.5, 1, 2];
us = linspace(1E-16,1-1E-16,1000);

for i = 1:length(s)
    gs = func_inf2uni(gauss,   us, s(i), 'logistic');
    ls = func_inf2uni(laplace, us, s(i), 'logistic');
    cs = func_inf2uni(cauchy,  us, s(i), 'logistic');
    
    subplot(length(s),3,(i-1)*3+1)
    plot(us, gs)
    if i == 1
        title('gauss, using logistic')
    end
    ylabel(['s=' num2str(s(i))])
    subplot(length(s),3,(i-1)*3+2)
    plot(us, ls)
    if i == 1
        title('laplace, using logistic')
    end
    subplot(length(s),3,(i-1)*3+3)
    plot(us, cs)
    if i == 1
        title('cauchy, using logistic')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 5;
f = figure(5);
set(f, 'position', [100, 100, 800, 800])
ys = linspace(-5, 5, 1000);

for i = 1:length(s)
    gs = func_inf2warp(gauss,   ys, a, s(i), 'logistic');
    ls = func_inf2warp(laplace, ys, a, s(i), 'logistic');
    cs = func_inf2warp(cauchy,  ys, a, s(i), 'logistic');
    
    subplot(length(s),3,(i-1)*3+1)
    plot(ys, gs)
    if i == 1
        title('gauss, using logistic')
    end
    ylabel(['s=' num2str(s(i))])
    subplot(length(s),3,(i-1)*3+2)
    plot(ys, ls)
    if i == 1
        title('laplace, using logistic')
    end
    subplot(length(s),3,(i-1)*3+3)
    plot(ys, cs)
    if i == 1
        title('cauchy, using logistic')
    end
end
