
sig = 0.3;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas

% Gaussian
ref   = @(u) erfinv(u*2-1)*sqrt(2);

poly1  = setup_oned(30, 'type', 'Fourier',  'domain', [-10,10]);
opt1   = ftt_options('method', 'AMEN', 'ng_flag', true, 'oned_ref', poly1, ...
    'err_tol', 1E-10, 'loc_err_tol', 1E-10, 'max_als', 20, ...
    'kick_rank', 3, 'max_rank', 50);


poly2  = setup_oned(30, 'type', 'Fourier',  'domain', [-4,4]);
opt2   = ftt_options('method', 'AMEN', 'ng_flag', true, 'oned_ref', poly2, ...
    'err_tol', 1E-10, 'loc_err_tol', 1E-10, 'max_als', 20, ...
    'kick_rank', 3, 'max_rank', 50);

beta = 10.^(-3:1/2:0);

ftts = {};
irts = {};

ftts{1} = build_ftt(@(u) ratio_fun_banana(u, dat, sig, 0, beta(1), []), 2, [], opt1);
irts{1} = build_irt(ftts{1});

for i = 2:length(beta)
    ftts{i} = build_ftt(@(u) ratio_fun_banana(u, dat, sig, beta(i-1), beta(i), irts), 2, [], opt2);
    irts{i} = build_irt(ftts{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
for i = 2:length(beta)
    
    rf = ratio_fun_banana(uts, dat, sig, beta(i-1), beta(i), irts(1:i-1));
    bf = fun_banana(xts, dat, sig, beta(i));
    
    h = figure(i);
    set(h, 'position', [100, 100, 800, 400])
    subplot(1,2,1)
    contour(xs, ys, reshape(bf(:), 100, 100))
    %axis([-2, 3, -5, 5])
    subplot(1,2,2)
    contour(us, us, reshape(rf(:), 100, 100))
    
end
%}

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init

n  = 100;
xs = linspace(-20, 20, n);
ys = linspace(-20, 40, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';


us = linspace(-2.2, 2.2, n);
[u1,u2] = meshgrid(us, us);
uts = [u1(:), u2(:)]';


bf = fun_banana(xts, dat, sig, beta(1));
rf = exp(-0.5*sum(uts.^2, 1));

h = figure(1);
set(h, 'position', [100, 100, 400, 400])
%contour(xs, ys, reshape(bf(:), n, n), 5, 'color', red, 'linewidth', 2)
contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\nu_0$', 'interpreter', 'latex', 'fontsize', 20)
colormap default


h = figure(2);
set(h, 'position', [100, 100, 400, 400])
contour(us, us, reshape(rf(:), n, n), 5, 'color', blue, 'linewidth', 2)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\mu$', 'interpreter', 'latex', 'fontsize', 20)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3
n  = 100;
xs = linspace(-5, 5, n);
ys = linspace(-5, 15, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';


us = linspace(-2.2, 2.2, n);
[u1,u2] = meshgrid(us, us);
uts = [u1(:), u2(:)]';

i  = 3;
bf = fun_banana(xts, dat, sig, beta(i));
rf = ratio_fun_banana(uts, dat, sig, beta(i-1), beta(i), irts(1:i-1));

h = figure(3);
set(h, 'position', [100, 100, 400, 400])
%contour(xs, ys, reshape(bf(:), n, n), 5, 'color', red, 'linewidth', 1)
contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\nu_1$', 'interpreter', 'latex', 'fontsize', 20)
colormap default


h = figure(4);
set(h, 'position', [100, 100, 400, 400])
%contour(us, us, reshape(rf(:), n, n), 5, 'color', purple, 'linewidth', 1)
contour(us, us, reshape(rf(:), n, n), 5, 'linewidth', 2)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('\vphantom{$\nu_0$}', 'interpreter', 'latex', 'fontsize', 20)
colormap summer



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5
n  = 100;
xs = linspace(-3, 3, n);
ys = linspace(-2, 6, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';


us = linspace(-2.2, 2.2, n);
[u1,u2] = meshgrid(us, us);
uts = [u1(:), u2(:)]';

i  = 5;
bf = fun_banana(xts, dat, sig, beta(i));
rf = ratio_fun_banana(uts, dat, sig, beta(i-1), beta(i), irts(1:i-1));

h = figure(5);
set(h, 'position', [100, 100, 400, 400])
%contour(xs, ys, reshape(bf(:), n, n), 5, 'color', red, 'linewidth', 1)
contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\nu_2$', 'interpreter', 'latex', 'fontsize', 20)
colormap default


h = figure(6);
set(h, 'position', [100, 100, 400, 400])
%contour(us, us, reshape(rf(:), n, n), 5, 'color', purple, 'linewidth', 1)
contour(us, us, reshape(rf(:), n, n), 5, 'linewidth', 2)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('\vphantom{$\nu_0$}', 'interpreter', 'latex', 'fontsize', 20)
colormap summer



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7
n  = 100;
xs = linspace(-1.5, 1.5, n);
ys = linspace(-1.2, 2, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';


us = linspace(-2.2, 2.2, n);
[u1,u2] = meshgrid(us, us);
uts = [u1(:), u2(:)]';

i  = 7;
bf = fun_banana(xts, dat, sig, beta(i));
rf = ratio_fun_banana(uts, dat, sig, beta(i-1), beta(i), irts(1:i-1));

h = figure(7);
set(h, 'position', [100, 100, 400, 400])
%contour(xs, ys, reshape(bf(:), n, n), 5, 'color', red, 'linewidth', 1)
contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\nu_3\equiv\nu_\pi$', 'interpreter', 'latex', 'fontsize', 20)
colormap default


h = figure(8);
set(h, 'position', [100, 100, 400, 400])
%contour(us, us, reshape(rf(:), n, n), 5, 'color', purple, 'linewidth', 1)
contour(us, us, reshape(rf(:), n, n), 5, 'linewidth', 2)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('\vphantom{$\nu_0$}', 'interpreter', 'latex', 'fontsize', 20)
colormap summer


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4
n  = 100;
xs = linspace(-4, 4, n);
ys = linspace(-3, 10, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';


us = linspace(-2.2, 2.2, n);
[u1,u2] = meshgrid(us, us);
uts = [u1(:), u2(:)]';

i  = 4;
bf = fun_banana(xts, dat, sig, beta(i));
rf = ratio_fun_banana(uts, dat, sig, beta(i-1), beta(i), irts(1:i-1));

h = figure(9);
set(h, 'position', [100, 100, 400, 400])
contour(xs, ys, reshape(bf(:), n, n), 5, 'color', red, 'linewidth', 1)
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\nu_1$', 'interpreter', 'latex', 'fontsize', 20)


h = figure(10);
set(h, 'position', [100, 100, 400, 400])
contour(us, us, reshape(rf(:), n, n), 5, 'color', purple, 'linewidth', 1)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('\vphantom{$\nu_0$}', 'interpreter', 'latex', 'fontsize', 20)

