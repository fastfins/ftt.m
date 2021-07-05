
% setup the OU process
d = 20;
a = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = setup_ou_process(d, a);

func1 = @(x) eval_ou_process(data, x);

%%%%

debug_size = 1E4;
debug_x = data.B\randn(d, debug_size);
sample_x = data.B\randn(d, 1E3);

%%%%

% setup the reference polynomial
poly1 = setup_oned(21, 'type', 'Chebyshev1st', 'domain', [-5,5]);
poly2 = setup_oned(21, 'type', 'Legendre', 'domain', [-5,5]);
poly3 = setup_oned(10, 'type', 'Fourier',  'domain', [-5,5]);
poly4 = setup_oned(5, 'type', 'Lagrange', 'lag_elems', 4, 'ghost_size', 1E-5, 'domain', [-5,5]);

opt1 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly1, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt1 = build_ftt(func1, d, [], opt1, 'debug_x', debug_x, 'sample_x', sample_x);

% squared version using AMEN
% AMEN has higher rank than random sampling with the same accuracy
opt2 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly2, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt2 = build_ftt(func1, d, [], opt2, 'debug_x', debug_x, 'sample_x', sample_x);
   

opt3 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly3, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt3 = build_ftt(func1, d, [], opt3, 'debug_x', debug_x, 'sample_x', sample_x);


opt4 = ftt_options('method', 'Random', 'ng_flag', true, 'oned_ref', poly4, ...
    'err_tol', 1E-4, 'loc_err_tol', 1E-10, 'max_rank', 19, 'max_als', 20);
ftt4 = build_ftt(func1, d, [], opt4, 'debug_x', debug_x, 'sample_x', sample_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = rand(d, 1E2);
debug_jac(ftt3, z, -1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
tol = 1E-3;
poly = poly1;
f   = rand(poly.num_nodes,1);
xs  = linspace(poly.domain(1)+tol*10, poly.domain(2)-tol*10, 100);
xsp = xs + tol;
fx  = eval_oned_nodes2x(poly, xs, f);
fxp = eval_oned_nodes2x(poly, xsp, f);
plot((fxp-fx)/tol, eval_oned_nodes2x_deri(poly, xs, f), '.');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function debug_jac(ftt, z, di)
ftt_int = build_irt(ftt, 'dir', di);
%
tic;ur = eval_irt(ftt_int, z);toc
%
tol = 1E-5;
tic;
[m,n] = size(ur);
ut = reshape(repmat(ur, m, 1), m, []) + repmat(eye(m)*tol, 1, n);
zt = eval_rt(ftt_int, ut);
%
um = reshape(repmat(ur, m, 1), m, []) - repmat(eye(m)*tol, 1, n);
zm = eval_rt(ftt_int, um);
%
Jd = ( zt - zm ) / (tol*2);
toc
%{
tic;
[m,n] = size(ur);
ut = reshape(repmat(ur, m, 1), m, []) + repmat(eye(m)*tol, 1, n);
zt = eval_rt(ftt_int, ut);
Jd = ( zt - reshape(repmat(z, m, 1), m, []) ) / tol;
toc
%}
%
tic; J = eval_rt_jac(ftt_int, ur, z); toc
%
norm(J(:)-Jd(:))
end