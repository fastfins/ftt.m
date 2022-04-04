
sig = 0.4;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas
fun = @(z) fun_banana(z, dat, sig, 1);
fun2 = @(z) fun_banana2(z, dat, sig, 1);

% Gaussian
%ref   = @(u) erfinv(u*2-1)*sqrt(2);

red    = '#D95319';
blue   = '#0072BD';
purple = '#7E2F8E';
grey = [0.7, 0.7, 0.7];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
nsteps = 2^10;
init = [1.2; 1.5];
tic;
out1 = NUTS(@(x) log_target(x,dat,sig,1), init, nsteps);
xx1 = out1.samples;
toc

init = [-0.9; -0.5];
tic;
out2 = NUTS(@(x) log_target(x,dat,sig,1), init, nsteps);
xx2 = out2.samples;
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diag = GaussReference();
poly = {Lagrangep(2,25,[-4,4]), Fourier(25, diag.domain)};
airt = DIRT(fun,2,poly,diag,'betas',[0.01, 0.15, 1]);

r = random(airt.ref, 2, 2E2);
xr = eval_irt(airt, r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n  = 10;
xxs = linspace(-2.5, 2.5, n);
xys = linspace(-2.5, 2.5, n);
[xx,yy] = meshgrid(xxs, xys);
s = [xx(:), yy(:)]';
%s = invert_cdf(airt.ref, s);
xs = eval_irt(airt, s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n  = 200;
rxs = linspace(airt.ref.domain(1), airt.ref.domain(2), n);
rys = linspace(airt.ref.domain(1), airt.ref.domain(2), n);
[xx,yy] = meshgrid(rxs, rys);
rts = [xx(:), yy(:)]';
%
ref = exp(-0.5*sum(rts.^2,1));
%
[mllkd, mlp] = fun(rts);
bf = exp(-mllkd-mlp);
const = 64/n^2;
bf = bf/(sum(bf(:))*const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position', [100, 100, 400, 400])

figure('position', [100, 100, 400, 400])
contour(rxs, rys, reshape(ref(:), n, n), 5, 'linewidth', 2, 'Color', grey)
hold on
plot(r(1,:),r(2,:), '.', 'Color', blue, 'markersize', 20);
%plot(s(1,:),s(2,:), '.', 'Color', red, 'markersize', 20);
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
%title('$\mu$', 'interpreter', 'latex', 'fontsize', 20)
axis([-2, 2, -2, 2])

figure('position', [100, 100, 400, 400])
contour(rxs, rys, reshape(bf(:), n, n), 5, 'linewidth', 2, 'Color', grey)
hold on
plot(xr(1,:),xr(2,:), '.', 'Color', blue, 'markersize', 20);
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
%title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
axis([-1.5, 1.5, -1, 2])

figure('position', [100, 100, 400, 400])
contour(rxs, rys, reshape(ref(:), n, n), 5, 'linewidth', 2, 'Color', grey)
hold on
plot(s(1,:),s(2,:), '.', 'Color', red, 'markersize', 20);
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
%title('$\mu$', 'interpreter', 'latex', 'fontsize', 20)
axis([-2, 2, -2, 2])

figure('position', [100, 100, 400, 400])
contour(rxs, rys, reshape(bf(:), n, n), 5, 'linewidth', 2, 'Color', grey)
hold on
plot(xs(1,:),xs(2,:), '.', 'Color', red, 'markersize', 20);
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
%title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
axis([-1.5, 1.5, -1, 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position', [100, 100, 400, 400])
contour(rxs, rys, reshape(bf(:), n, n), 5, 'linewidth', 2, 'Color', grey)
hold on
plot(xx1(1,1:50),xx1(2,1:50), '.-', 'Color', red,  'markersize', 10, 'LineWidth', 0.5);
plot(xx2(1,1:50),xx2(2,1:50), '.-', 'Color', blue, 'markersize', 10, 'LineWidth', 0.5);
set(gca, 'fontsize', 20, 'TickLabelInterpreter','latex')
%title('$\pi$', 'interpreter', 'latex', 'fontsize', 20)
axis([-1.5, 1.5, -1, 2])

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

function [mlf,gmlf] = log_target(u,data,sigma,beta)

F   = log( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
mllkd = sum((F-data).^2,1)*beta/(2*sigma^2);
mlp = 0.5*sum(u.^2,1);

mlf = mllkd + mlp;

gmlp = u;
tmp = [2*(u(1,:)-1)+400*(u(1,:).^2-u(2,:)).*u(1,:); 200*(u(2,:)-u(1,:).^2)]...
    ./( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
gmllkd = (F-data(1)).*tmp*beta/sigma^2 + (F-data(2)).*tmp*beta/sigma^2;

gmlf = gmllkd+gmlp;

end
