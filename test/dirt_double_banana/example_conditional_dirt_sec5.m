close all

% define the joint density
sig = 0.3;
fun = @(z) joint_banana(z,sig);

refmap = GaussMap([-4, 4]);
poly = Fourier(25, [-4, 4]); 
opt = FTToption('max_als', 2, 'als_tol', 1E-8, 'local_tol', 1E-5, 'kick_rank', 2, 'init_rank', 20, 'max_rank', 30);
if ~exist('irt')
    irt = DIRT(fun, 3, poly, refmap, opt, 'ess_tol', 0.5); % the conditonal DIRT
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data = [-2, -1, 0, 1, 2, 3];

data = 0;

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';
grey = [0.7, 0.7, 0.7];

n  = 200;
xs = linspace(-4, 4, n);
ys = linspace(-4, 4, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';
%
us = linspace(0, 1, n);
rxs = eval_icdf(irt.diag, us);
rys = rxs;
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';

for ii = 1:length(data)
    dat = data(ii); % data
    ry = eval_rt(irt, dat); % reference data sample
    
    % true conditional, unnormalised
    [mllkd,mlp] = fun([repmat(dat,1,size(xts,2));xts]);
    bf = exp(-mllkd-mlp);
    % conditional irt density in target space, unnormalised
    rf = log_pdf(irt, [repmat(dat,1,size(xts,2));xts]);
    
    figure('position', [100, 100, 1200, 800])
    subplot(2,4,1)
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    subplot(2,4,5)
    contour(xs, ys, reshape(exp(rf(:)), n, n), 10, 'linewidth', 1)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    % pullback density of the true conditinal
    mlf = pullback(irt, fun, [repmat(ry,1,size(rts,2));rts]);
    %
    
    r = random(irt.diag, 2, 256); % sample the reference measure
    subplot(2,4,2)
    %contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 20, 'linewidth', 1)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 10, 'linewidth', 1, 'Color', grey)
    hold on
    plot(r(1,:), r(2,:), '.', 'MarkerSize', 15, 'Color', purple)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$T^\sharp \pi$, random', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    
    xx = eval_irt(irt, [repmat(ry,1,size(r,2));r]);
    xx = xx(2:3,:);
    subplot(2,4,6)
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', 1, 'Color', grey)
    hold on
    plot(xx(1,:), xx(2,:), '.', 'MarkerSize', 15, 'Color', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    
    r = sobol(irt.diag, 2, 256); % QMC sample the reference measure
    subplot(2,4,3)
    %contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 20, 'linewidth', 1)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 10, 'linewidth', 1, 'Color', grey)
    hold on
    plot(r(1,:), r(2,:), '.', 'MarkerSize', 15, 'Color', purple)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$T^\sharp \pi$, Sobol', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    
    xx = eval_irt(irt, [repmat(ry,1,size(r,2));r]);
    xx = xx(2:3,:);
    subplot(2,4,7)
    contour(xs, ys, reshape(exp(rf(:)), n, n), 10, 'linewidth', 1, 'Color', grey)
    hold on
    plot(xx(1,:), xx(2,:), '.', 'MarkerSize', 15, 'Color', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    % regular grid in the reference measure
    [xx,yy] = meshgrid(linspace(-4,4,16), linspace(-4,4,16));
    r = [xx(:), yy(:)]';
    subplot(2,4,4)
    %contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 20, 'linewidth', 1)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 10, 'linewidth', 1, 'Color', grey)
    hold on
    plot(r(1,:), r(2,:), '.', 'MarkerSize', 15, 'Color', purple)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$T^\sharp \pi$, grid', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    
    xx = eval_irt(irt, [repmat(ry,1,size(r,2));r]);
    xx = xx(2:3,:);
    subplot(2,4,8)
    contour(xs, ys, reshape(exp(rf(:)), n, n), 10, 'linewidth', 1, 'Color', grey)
    hold on
    plot(xx(1,:), xx(2,:), '.', 'MarkerSize', 15, 'Color', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    
    nsteps = 2^10;
    init = [0.9; 0];
    tic;
    out1 = NUTS(@(x) log_target(fun,dat,x), init, nsteps);
    xx1 = out1.samples;
    toc
    %
    tic;
    out2 = NUTS(@(z) log_target_pullback_nuts(irt,fun,ry,z), init, nsteps);
    xs2 = eval_irt(irt, [repmat(ry,1,size(out2.samples,2));out2.samples]);
    xx2 = xs2(2:3,:);
    toc
    tic; % Metropolized independence sampler
    out3 = pCN(@(z) log_target_pullback_pcn(irt,fun,ry,z), init, nsteps, log(10));
    xs3 = eval_irt(irt, [repmat(ry,1,size(out3.samples,2));out3.samples]);
    xx3 = xs3(2:3,:);
    toc
    tic; % negatively correlated pCN
    out4 = pCN(@(z) log_target_pullback_pcn(irt,fun,ry,z), init, nsteps, log(2));
    xs4 = eval_irt(irt, [repmat(ry,1,size(out4.samples,2));out4.samples]);
    xx4 = xs4(2:3,:);
    toc   
    
    msteps = 1:5;
    figure('position', [100, 100, 1200, 600])
    subplot(4,4,[1,5])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', [.8 .8 .8])
    hold on
    plot(xx1(1,:), xx1(2,:), '.', 'Color', blue)
    plot(xx1(1,msteps), xx1(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('NUTS', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
        
    subplot(4,4,[2,6])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', [.8 .8 .8])
    hold on
    plot(xx2(1,:), xx2(2,:), '.', 'Color', blue)
    plot(xx2(1,msteps), xx2(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('preconditioned NUTS', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    subplot(4,4,[3,7])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', [.8 .8 .8])
    hold on
    plot(xx3(1,:), xx3(2,:), '.', 'Color', blue)
    plot(xx3(1,msteps), xx3(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('Metropilis independence', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    
    subplot(4,4,[4,8])
    contour(xs, ys, reshape(bf(:), n, n), 10, 'linewidth', .5, 'Color', [.8 .8 .8])
    hold on
    plot(xx4(1,:), xx4(2,:), '.', 'Color', blue)
    plot(xx4(1,msteps), xx4(2,msteps), '-o', 'linewidth', 1, 'Color', purple, 'MarkerFaceColor', purple)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 16)
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    title('preconditioned pCN', 'interpreter', 'latex', 'fontsize', 16)
    colormap default
    axis([-2, 2, -1, 2])
    
    subplot(4,4,9);  autocorr(xx1(1,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel(''); title('')
    subplot(4,4,13); autocorr(xx2(2,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    %     
    subplot(4,4,10); autocorr(xx2(1,:));  
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel(''); title('')
    subplot(4,4,14); autocorr(xx2(2,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    %     
    subplot(4,4,11); autocorr(xx3(1,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel(''); title('')
    subplot(4,4,15); autocorr(xx3(2,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    %
    subplot(4,4,12); autocorr(xx4(1,:)); 
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_1$','interpreter', 'latex', 'fontsize', 16); 
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    xlabel(''); title('')
    subplot(4,4,16); autocorr(xx4(2,:));
    set(gca, 'fontsize', 16, 'TickLabelInterpreter','latex')
    ylabel('$x_2$','interpreter', 'latex', 'fontsize', 16)
    xlabel('Lag','interpreter', 'latex', 'fontsize', 16)
    title('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mlf,gmlf] = log_target(func,y,x)

[mllkd,mlp,gmllkd,gmlp] = func([repmat(y,1,size(x,2));x]);
%
mlf = mllkd + mlp;
ind = length(y) + (1:size(x,1));
gmlf = gmllkd(ind,:)+gmlp(ind,:);

end

function [mlf,gmlf] = log_target_pullback_nuts(irt,func,ry,z)

[mlf,gmlf] = pullback(irt, func, [repmat(ry,1,size(z,2));z]);
%
ind = length(ry) + (1:size(z,1));
gmlf = gmlf(ind,:);

end

function [mllkd,mlp] = log_target_pullback_pcn(irt,func,ry,z)

mlf = pullback(irt, func, [repmat(ry,1,size(z,2));z]);
%
mlp = 0.5*sum(z.^2,1);
mllkd = mlf - mlp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f1, f2, g1, g2] = joint_banana(z, sigma)

y = z(1,:);
u = z(2:3,:);
F = log( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
G = [F-3; F-5];
%
mllkd = sum((G-y).^2,1)/(2*sigma^2);
mlp = 0.5*sum(u.^2,1);

f1 = mllkd + mlp - 0.5*sum(z.^2,1);
f2 = 0.5*sum(z.^2,1);

if nargout > 2
    tmp = [2*(u(1,:)-1)+400*(u(1,:).^2-u(2,:)).*u(1,:); 200*(u(2,:)-u(1,:).^2)]...
        ./( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
    g1 = [-(2*F-2*y-8); (2*F-2*y-8)*tmp]/sigma^2;
    g2 = z;
end

end