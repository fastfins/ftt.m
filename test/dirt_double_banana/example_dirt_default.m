
sig = 0.2;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas
fun = @(z) fun_banana(z, dat, sig, 1);

% Gaussian
ref   = @(u) erfinv(u*2-1)*sqrt(2);

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';

% betas = 0.01, 0.125, 0.8, 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airt = DIRT(fun,2,[-4,4]);

n  = 100;
rxs = linspace(airt.ref.domain(1), airt.ref.domain(2), n);
rys = linspace(airt.ref.domain(1), airt.ref.domain(2), n);
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
        logfz = log_joint_pdf(airt.ref, rts);
        
        bf = exp(-mllkd*(airt.betas(k)-airt.betas(k-1))+logfz);
        
        subplot(2,2,4)
        contour(rxs, rys, reshape(bf(:), n, n), 8, 'linewidth', 1)
        xlabel('$u_1$', 'interpreter', 'latex', 'fontsize', 20)
        ylabel('$u_2$', 'interpreter', 'latex', 'fontsize', 20)
        set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
        title('a. ratio fun', 'interpreter', 'latex', 'fontsize', 20)
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
