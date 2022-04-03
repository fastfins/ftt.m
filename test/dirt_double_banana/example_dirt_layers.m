
sig = 0.2;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas
fun = @(z) fun_banana(z, dat, sig, 1);
fun2 = @(z) fun_banana2(z, dat, sig, 1);

% Gaussian
%ref   = @(u) erfinv(u*2-1)*sqrt(2);

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';

% betas = 0.01, 0.125, 0.8, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target

%n  = 100;
%xs = linspace(-4, 4, n);
%ys = linspace(-4, 4, n);
%[xx,yy] = meshgrid(xs, ys);
%xts = [xx(:), yy(:)]';
%const = 64/n^2;

%[mllkd, mlp] = fun(xts);
%bf = exp(-mllkd-mlp);
%rf = exp(-0.5*sum(xts.^2, 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diag = GaussReference();
poly = {Lagrangep(2,25,[-4,4]), Fourier(25, diag.domain)};
airt = DIRT(fun,2,poly,diag,'betas',[0.004, 0.05, 0.4, 1]);

debug_r = random(airt, 1E4);
sample_r = random(airt, 1E4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = FTToption('tt_method', 'amen', 'sqrt_flag', true, ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 50, 'max_als', 10);
irt = SIRT(fun2, 2, Lagrange1(1000,[-4,4]), opt, 'sample_x', sample_r);

%f = eval_pdf(irt, debug_r);
irt_r = random(irt, 5E3);
f = eval_pdf(irt, debug_r);
f2 = fun2(debug_r);

% squared Hellinger error between logfx and (-mllkds-mlps)
[~,dh2,~] = f_divergence(log(f), log(f2));
disp(['Hellinger error of single layer (1000 points per dimension): ' num2str(sqrt(dh2))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n  = 400;
rxs = linspace(airt.ref.domain(1), airt.ref.domain(2), n);
rys = linspace(airt.ref.domain(1), airt.ref.domain(2), n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';

figure('position', [100, 100, 400, 400])

figure('position', [100, 100, 400, 400])
ref = exp(-0.5*sum(rts.^2,1));
contour(rxs, rys, reshape(ref(:), n, n), 5, 'linewidth', 2, 'Color', blue)
%xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
%ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\mu$', 'interpreter', 'latex', 'fontsize', 20)
axis([-2, 2, -2, 2]*1.1)

for k = 1:airt.n_layers
    figure('position', [100, 100, 400, 400])
    [mllkd, mlp] = fun(rts);
    bf = exp(-mllkd*airt.betas(k)-mlp);
    const = 64/n^2;
    bf = bf/(sum(bf(:))*const);
    contour(rxs, rys, reshape(bf(:), n, n), 5, 'linewidth', 2)
    %xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    %ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    if k == airt.n_layers
    title(['$\pi_' num2str(k-1) '=\pi$'], 'interpreter', 'latex', 'fontsize', 20)
    else
    title(['$\pi_' num2str(k-1) '$'], 'interpreter', 'latex', 'fontsize', 20)
    end
    colormap default
    axis([-2, 2, -2, 2]*1.1)
end

for k = 2:airt.n_layers
    [x,logf] = eval_irt(airt, rts, k-1);
    [mllkd, mlp] = fun(x);
    logfz = log_joint_pdf(airt.ref, rts);

    bf = exp(-mllkd*(airt.betas(k)-airt.betas(k-1))+logfz);

    figure('position', [100, 100, 400, 400])
    contour(rxs, rys, reshape(bf(:), n, n), 10, 'linewidth', 1)
    %xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    %ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    axis([-2, 2, -2, 2]*1.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position', [100, 100, 400, 400])
ref = exp(-0.5*sum(rts.^2,1));
contour(rxs, rys, reshape(ref(:), n, n), 5, 'linewidth', 1, 'Color', blue)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\mu$', 'interpreter', 'latex', 'fontsize', 20)
axis([-2, 2, -2, 2]*1.1)

figure('position', [100, 100, 400, 400])
[mllkd, mlp] = fun(rts);
bf = exp(-mllkd-mlp);
const = 64/n^2;
bf = bf/(sum(bf(:))*const);
contour(rxs, rys, reshape(bf(:), n, n), 2, 'linewidth', 1, 'Color', blue)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
axis([-2, 2, -2, 2]*1.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mllkd, mlp] = fun_banana(u, data, sigma, beta)

F   = log((1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2);
mllkd = sum((F-data).^2,1)*beta/(2*sigma^2);
mlp = 0.5*sum(u.^2,1);

%lpt = -sum((F-data).^2,1)*beta/(2*sigma^2) - sum(u.^2,1)*beta/2;%
%p   = exp(lpt);

end


function f = fun_banana2(u, data, sigma, beta)

F   = log((1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2);
mllkd = sum((F-data).^2,1)*beta/(2*sigma^2);
mlp = 0.5*sum(u.^2,1);
f = exp(-mllkd-mlp);
%lpt = -sum((F-data).^2,1)*beta/(2*sigma^2) - sum(u.^2,1)*beta/2;%
%p   = exp(lpt);

end
