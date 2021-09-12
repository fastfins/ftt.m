
% setup the reference
tol = 0;
poly1 = Legendre(20, [tol, 1]);
poly2 = Lagrange1(50);
poly3 = Lagrangep(5, 4, [tol, 1], 'ghost_size', 1E-10);

%%%%

d = 10;

func = @(x) sqrt(1./sum(1E-5+x.^2,1));

%%%%


debug_size = 1E4;
debug_x = zeros(d, debug_size);
for k = 1:d
    debug_x(k,:) = sample_domain(poly1, debug_size);
end

opt1 = FTToption('max_als', 4, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
opt2 = FTToption('max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);

polys = {poly1, poly2, poly3};
for i = 1:3
    tic;
    ftt{i,1} = FTT(func, d, polys{i}, opt1, 'debug_x', debug_x);
    ftt{i,2} = round(ftt{i,1}, 1E-4);
    
    ftt{i,3} = FTT(func, d, polys{i}, opt2, 'debug_x', debug_x);
    ftt{i,4} = round(ftt{i,3}, 1E-4);
    toc
end

x = debug_x(:,10);
for i = 1:3
    for j= 1:4
        [gx,fx] = grad(ftt{i,j}, x);
        f = eval(ftt{i,j}, x);
        tol = 1E-5;
        fp = zeros(size(gx));
        fm = zeros(size(gx));
        for ii = 1:length(x)
            xp = x;
            xp(ii) = xp(ii)+tol;
            xm = x;
            xm(ii) = xm(ii)-tol;
            fp(ii) = eval(ftt{i,j}, xp);
            fm(ii) = eval(ftt{i,j}, xm);
        end
        disp([norm(f-fx) norm((fp-fm)/(2*tol) - gx)])
    end
end


