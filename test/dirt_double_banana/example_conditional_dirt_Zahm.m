close all

% define the joint density
sig = 0.3;
fun = @(z) joint_banana(z,sig);

refmap = GaussMap([-4, 4]);
poly = Fourier(20, [-4, 4]); 
opt = FTToption('max_als', 2, 'als_tol', 1E-8, 'local_tol', 1E-5, 'kick_rank', 1, 'init_rank', 20, 'max_rank', 20);
if ~exist('irt')
    irt = DIRT(fun, 3, poly, refmap, opt, 'ess_tol', 0.5); % the conditonal DIRT
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data = [-2, -1, 0, 1, 2, 3];

data = 0;

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';

n  = 100;
xs = linspace(-4, 4, n);
ys = linspace(-4, 4, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';
%
rxs = linspace(refmap.domain(1), refmap.domain(2), n);
rys = linspace(refmap.domain(1), refmap.domain(2), n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';

for ii = 1:length(data)
    dat = data(ii); % data
    ry = eval_rt(irt, dat); % reference sample
    
    % true conditional, unnormalised
    [mllkd,mlp] = fun([repmat(dat,1,size(xts,2));xts]);
    bf = exp(-mllkd-mlp);
    % conditional irt density in target space, unnormalised
    rf = log_pdf(irt, [repmat(dat,1,size(xts,2));xts]);
    
    figure('position', [100, 100, 800, 800])
    subplot(2,2,1)
    contour(xs, ys, reshape(bf(:), n, n), 5, 'linewidth', 2)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    subplot(2,2,2)
    contour(xs, ys, reshape(exp(rf(:)), n, n), 5, 'linewidth', 2)
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    % conditional irt density in reference space, unnormalised
    [x,logf] = eval_irt(irt, [repmat(ry,1,size(rts,2));rts]);
    % pullback density of the true conditinal
    mlf = pullback(irt, fun, [repmat(ry,1,size(rts,2));rts]);
    %
    subplot(2,2,3)
    contour(rxs, rys, reshape(exp(logf(:)), n, n), 5, 'linewidth', 2)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$\hat\pi(u)$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    subplot(2,2,4)
    contour(rxs, rys, reshape(exp(-mlf(:)), n, n), 5, 'linewidth', 2)
    xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('$T^\sharp \hat\pi$', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    
    init = randn(2,1);
    tic;
    for irun=1:16
        out1 = NUTS(@(x) log_target(fun,dat,x), init, 2^12);
        xx1(:,:,irun) = out1.samples;
    end
    toc
    %
    tic;
    for irun=1:16
        out2 = NUTS(@(z) log_target_pullback_nuts(irt,fun,ry,z), init, 2^12);
        xs2 = eval_irt(irt, [repmat(ry,1,size(out2.samples,2));out2.samples]);
        xx2(:,:,irun) = xs2(2:3,:);
    end
    toc
    tic;
    for irun=1:16
        out3 = pCN(@(z) log_target_pullback_pcn(irt,fun,ry,z), init, 2^12, log(10));
        xs3 = eval_irt(irt, [repmat(ry,1,size(out3.samples,2));out3.samples]);
        xx3(:,:,irun) = xs3(2:3,:);
    end
    toc
    tic;
    for irun=1:16
        out4 = pCN(@(z) log_target_pullback_pcn(irt,fun,ry,z), init, 2^12, log(2));
        xs4 = eval_irt(irt, [repmat(ry,1,size(out4.samples,2));out4.samples]);
        xx4(:,:,irun) = xs4(2:3,:);
    end
    toc    
    
    figure('position', [100, 100, 1200, 800])
    subplot(4,3,[1,4])
    plot(out1.samples(1,:), out1.samples(2,:), '.')
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('HMC', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    axis([-4, 4, -4, 4])
        
    subplot(4,3,[2,5])
    plot(xx2(1,:), xx2(2,:), '.')
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('HMC, preconditioned', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    axis([-4, 4, -4, 4])
    
    subplot(4,3,[3,6])
    plot(xx3(1,:), xx3(2,:), '.')
    xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 20)
    set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
    title('pCN, preconditioned', 'interpreter', 'latex', 'fontsize', 20)
    colormap default
    axis([-4, 4, -4, 4])
    
%     subplot(4,3,7);  autocorr(out1.samples(1,:)); title('$x_1$','interpreter', 'latex'); 
%     subplot(4,3,10); autocorr(out1.samples(2,:)); title('$x_2$','interpreter', 'latex')
%     
%     subplot(4,3,8);  autocorr(xx2(1,:)); title('$x_1$','interpreter', 'latex'); 
%     subplot(4,3,11); autocorr(xx2(2,:)); title('$x_2$','interpreter', 'latex')
%     
%     subplot(4,3,9);  autocorr(xx3(1,:)); title('$x_1$','interpreter', 'latex'); 
%     subplot(4,3,12); autocorr(xx3(2,:)); title('$x_2$','interpreter', 'latex')
    
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