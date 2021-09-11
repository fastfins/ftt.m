d = 5; a = 0.5;
A = diag(-sqrt(1-a^2)*ones(d-1,1), -1) + eye(d);
D = diag([1, a*ones(1,d-1)]);
B = D\A;
Q = B'*B;   % precision matrix
C = inv(Q); % covariance matrix
z = sqrt((2*pi)^d/det(Q)); % normalising constant
% The joint distribution, unnormalised


nsteps = 1E4;
tic;
out = NUTS(@(x) log_target(B,Q,x), randn(d,1), nsteps);
toc

norm(cov(out.samples') - C, 'fro')/norm(C, 'fro')
mean(out.samples')

%%%%%

function [f,g] = log_target(B,Q,x)

f = 0.5*sum((B*x).^2,1);
g = Q*x;

end