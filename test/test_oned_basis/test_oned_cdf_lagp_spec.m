 

pdf1 = @(x) exp(-0.5*x.^2);
pdf2 = @(x) exp(-0.5*abs(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

lag = Lagrangep(6, 4, [-4, 4]);
lag_cdf = LagrangepCDF(lag);

fl = pdf2(lag.nodes);
pl = eval(lag, fl, lag_cdf.nodes);

xs = linspace(lag.domain(1), lag.domain(2), 1000);
fi = eval(lag, fl, xs);
Fi = eval_cdf(lag_cdf, pl, xs);

z = rand(5E5,1);
tic;r = invert_cdf(lag_cdf, pl, z);toc
tic; norm(eval_cdf(lag_cdf, pl, r)-z), toc

dl = pdf2cdf(lag_cdf, pl);

figure
plot(xs, fi/dl.norm, xs, Fi)
hold on
plot(lag_cdf.nodes, dl.cdf_nodes/dl.norm)
ksdensity(r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lag = Lagrangep(6, 5, [-5, 5], 'ghost_size', 1);
lag_cdf = LagrangepCDF(lag);

fl = pdf2(lag.nodes).^0.5;
pl = eval(lag, fl, lag_cdf.nodes).^2;

xs = linspace(lag.domain(1), lag.domain(2), 10000);
fi = eval(lag, fl, xs).^2;
Fi = eval_cdf(lag_cdf, pl, xs);

z = rand(5E5,1);
tic;r = invert_cdf(lag_cdf, pl, z);toc
tic; norm(eval_cdf(lag_cdf, pl, r)-z), toc

dl = pdf2cdf(lag_cdf, pl);

figure
plot(xs, fi/dl.norm, xs, Fi)
hold on
plot(lag.nodes, fl.^2/dl.norm, 'o', lag_cdf.nodes, pl/dl.norm, 's')
ksdensity(r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


poly = Legendre(8, [-6, 6]);
poly_cdf = ChebyshevCDF(poly);

fp = poly.node2basis*pdf2(poly.nodes);
pp = eval(poly, fp, poly_cdf.nodes);

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), 10000);
fi = eval(poly, fp, xs);
Fi = eval_cdf(poly_cdf, pp, xs);

figure
plot(xs, fi, xs, Fi)
hold on
plot(poly_cdf.sampling_nodes, cdf_nodes, 'o')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = Legendre(8, [-6, 6]);
poly_cdf = ChebyshevCDF(poly);

fp = poly.node2basis*( pdf2(poly.nodes).^(0.5) );
pp = eval(poly, fp, poly_cdf.nodes).^2;
%
cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), 1000);
fi = eval(poly, fp, xs).^2;
Fi = eval_cdf(poly_cdf, pp, xs);


z = rand(5E5,1);
tic;r = invert_cdf(poly_cdf, pp, z);toc
tic; norm(eval_cdf(poly_cdf, pp, r)-z), toc

figure
plot(xs, fi, xs, Fi)
hold on
plot(poly_cdf.sampling_nodes, cdf_nodes, 'o')
ksdensity(r)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multi cdf test
n = 5E5;

lag = Lagrangep(6, 4, [-5, 5], 'ghost_size', 1);
lag_cdf = LagrangepCDF(lag);

fl = pdf2(lag.nodes).^(0.5);
pl = eval(lag, fl, lag_cdf.nodes).^2;
%
pl = repmat(pl, 1, n);
dl = pdf2cdf(lag_cdf, pl);

xs = linspace(lag.domain(1), lag.domain(2), n);
fi = eval(lag, fl, xs).^2;
Fi = eval_cdf(lag_cdf, pl, xs);

z = rand(n,1);
tic;r = invert_cdf(lag_cdf, pl, z);toc
tic; norm(eval_cdf(lag_cdf, pl, r)-z), toc

figure
plot(xs, fi/dl.norm(1), xs, Fi)
hold on
plot(lag.nodes, fl.^2/dl.norm(1), 'o', lag_cdf.nodes, pl(:,1)/dl.norm(1), 's')
ksdensity(r)


%%%%%%

poly = Legendre(8, [-6, 6]);
poly_cdf = ChebyshevCDF(poly);

fp = poly.node2basis*( pdf2(poly.nodes).^(0.5) );
pp = eval(poly, fp, poly_cdf.nodes).^2;
pp = repmat(pp, 1, n);

cdf_nodes = poly_cdf.cdf_basis2node*(poly_cdf.node2basis*pp);
cdf_nodes = cdf_nodes - cdf_nodes(1,:);
cdf_nodes = cdf_nodes./cdf_nodes(end,:);

xs = linspace(poly.domain(1), poly.domain(2), n);
fi = eval(poly, fp, xs).^2;
Fi = eval_cdf(poly_cdf, pp, xs);

z = rand(n,1);
tic;r = invert_cdf(poly_cdf, pp, z);toc
tic; norm(eval_cdf(poly_cdf, pp, r)-z), toc

figure
plot(xs, fi.^2, xs, Fi)
hold on
plot(poly_cdf.nodes, pp(:,1), 'o', poly_cdf.sampling_nodes, cdf_nodes(:,1), 's')
ksdensity(r)

