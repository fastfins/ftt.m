
% setup the reference
tol = 0;
poly1 = Legendre(20, [tol, 1]);
poly2 = Lagrange1(50);
poly3 = Lagrangep(5, 4, [tol, 1], 'ghost_size', 1E-10);

%%%%

d = 10;

func1 = @(x) sqrt(1./sum(1E-5+x.^2,1));
func2 = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
func3 = @(x) [(1 + sum(x,1)).^(-d-1); exp( - sum(abs(x - 0.5), 1)); cos( sum(x.^2,1) )];

func = func2;

%%%%


debug_size = 1E4;
debug_x = zeros(d, debug_size);
for k = 1:d
    debug_x(k,:) = sample_domain(poly1, debug_size);
end

opt1 = FTToption('max_als', 6, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
opt2 = FTToption('tt_method', 'random', 'max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-10, 'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);

polys = {poly1, poly2, poly3};
for i = 1:3
    tic;
    ftt{i,1} = FTT(func, d, polys{i}, opt1, 'debug_x', debug_x);
    ftt{i,2} = round(ftt{i,1}, 1E-4);
    ftt{i,3} = FTT(func, d, ftt{i,1}, 'debug_x', debug_x);
    ftt{i,4} = FTT(func, d, ftt{i,2}, 'debug_x', debug_x);
    
    ftt{i,5} = FTT(func, d, polys{i}, opt2, 'debug_x', debug_x);
    ftt{i,6} = round(ftt{i,5}, 1E-4);
    ftt{i,7} = FTT(func, d, ftt{i,5}, 'debug_x', debug_x);
    ftt{i,8} = FTT(func, d, ftt{i,6}, 'debug_x', debug_x);
    toc
end

%{
if FTT.direction > 0
    plot_ftt_1d(FTT, func, d)
else
    plot_ftt_1d(FTT, func, 1)
end
%}

%{
figure
[xx, yy] = meshgrid(linspace(0, 1, 100), linspace(0, 1, 100));
f = reshape(sqrt(1./(xx(:).^2 + yy(:).^2 + 0.001^2)), 100, 100);
surf(f)
%}
%
tic; exact   = func(debug_x); toc
for i = 1:3
    for j = 1:8
        tic; approx{i,j} = eval(ftt{i,j}, debug_x); toc
    end
end

for i = 1:3
    figure
    for j = 1:8
        plot(exact(:) - approx{i,j}(:), '.')
        hold on
    end
    legend('amen', 'amen rounded', 'amen rebuild', 'amen rounded rebuild', 'rand', 'rand rounded' , 'random rebuild', 'random rounded rebuild')
end

