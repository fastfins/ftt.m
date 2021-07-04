
function [J, Jd] = debug_jac(ftt, z, di)
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
%
tic; J = eval_rt_jac(ftt_int, ur, z); toc
%
norm(J(:)-Jd(:))

plot(J(:), Jd(:), '.')
end