
sig = 0.3;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas
fun = @(z) fun_banana(z, dat, sig, 1);

% Gaussian
ref   = @(u) erfinv(u*2-1)*sqrt(2);

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target

n  = 100;
xs = linspace(-4, 4, n);
ys = linspace(-4, 4, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';
const = 64/n^2;

[mllkd, mlp] = fun(xts);
bf = exp(-mllkd-mlp);
rf = exp(-0.5*sum(xts.^2, 1));

%{
h = figure(1);
set(h, 'position', [100, 100, 800, 400])
subplot(1,2,1)
contour(xs, ys, reshape(rf(:), n, n), 5, 'color', blue, 'linewidth', 2)
xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\mu$', 'interpreter', 'latex', 'fontsize', 20)
subplot(1,2,2)
contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
colormap default
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diag = UniformMap();
%diag = GaussMap([-4, 4]);

poly1 = {Legendre(60, [-4, 4]), Legendre(40, diag.domain)};
poly2 = {Lagrange1(60, [-4, 4], 'ghost_size', 1E-2, 'bc', 'Dirichlet'), Lagrange1(100, diag.domain, 'ghost_size', 1E-2, 'bc', 'Dirichlet')};
poly3 = {Lagrangep(5, 12, [-4, 4]), Lagrangep(5, 10, diag.domain, 'ghost_size', 1E-2, 'bc', 'Dirichlet')};
poly4 = {Fourier(30, [-4, 4]), Fourier(30, diag.domain)};

opt1 = FTToption('max_als', 1, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 40, 'max_rank', 50);
opt2 = FTToption('tt_method', 'random', 'max_als', 1, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 40, 'max_rank', 50);

%irt = DIRT(fun, 2, poly2, diag, opt1, 'min_beta', 1E-3, 'ess_tol', 0.8, 'betas',  2.^(-9:0));

airt = DIRT(fun, 2, poly3, diag, opt1, 'min_beta', 1E-3, 'ess_tol', 0.8, 'method', 'Aratio');
eirt = DIRT(fun, 2, poly3, diag, opt1, 'min_beta', 1E-3, 'ess_tol', 0.8, 'method', 'Eratio');


n  = 100;
rxs = linspace(diag.domain(1), diag.domain(2), n);
rys = linspace(diag.domain(1), diag.domain(2), n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';

for k = 1:airt.n_layers
    figure('position', [100, 100, 1200, 800])
    
    rf = log_pdf(airt, xts, k);
    rf = exp(rf);
    [mllkd, mlp] = fun(xts);
    bf = exp(-mllkd*airt.betas(k)-mlp);
    bf = bf/(sum(bf(:))*const);
    
    subplot(2,2,1)
    contour(xs, ys, reshape(rf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    
    subplot(2,2,2)
    contour(xs, ys, reshape(bf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    fk = eval_pdf(airt.irts{k}, rts);
    subplot(2,2,3)
    contour(rxs, rys, reshape(fk(:), n, n), 8, 'linewidth', 1)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('tt', 'interpreter', 'latex', 'fontsize', 20)
    
    if k > 1
        [x,logf] = eval_irt(airt, rts, k-1);
        [mllkd, mlp] = fun(x);
        logfz = log_pdf(diag, rts);
        
        bf = exp(-mllkd*(airt.betas(k)-airt.betas(k-1))+logfz);
        
        subplot(2,2,4)
        contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
        xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
        ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
        set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
        title('a. ratio fun', 'interpreter', 'latex', 'fontsize', 20)
    end
end


for k = 1:eirt.n_layers
    figure('position', [100, 100, 1200, 800])
    
    rf = log_pdf(eirt, xts, k);
    rf = exp(rf);
    [mllkd, mlp] = fun(xts);
    bf = exp(-mllkd*eirt.betas(k)-mlp);
    bf = bf/(sum(bf(:))*const);
    
    subplot(2,2,1)
    contour(xs, ys, reshape(rf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    
    subplot(2,2,2)
    contour(xs, ys, reshape(bf(:), n, n), 8, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    fk = eval_pdf(eirt.irts{k}, rts);
    subplot(2,2,3)
    contour(rxs, rys, reshape(fk(:), n, n), 8, 'linewidth', 1)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('tt', 'interpreter', 'latex', 'fontsize', 20)
    
    if k > 1
        [x,logf] = eval_irt(eirt, rts, k-1);
        [mllkd, mlp] = fun(x);
        logfz = log_pdf(diag, rts);
        
        bf = exp(-mllkd*eirt.betas(k)-mlp-logf+logfz);
        
        subplot(2,2,4)
        contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
        xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
        ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
        set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
        title('e. ratio fun', 'interpreter', 'latex', 'fontsize', 20)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mllkd, mlp] = fun_banana(u, data, sigma, beta)

F   = log((1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2);
mllkd = sum((F-data).^2,1)*beta/(2*sigma^2);
mlp = 0.5*sum(u.^2,1);

%lpt = -sum((F-data).^2,1)*beta/(2*sigma^2) - sum(u.^2,1)*beta/2;%
%p   = exp(lpt);

end
