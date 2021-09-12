function out = pCN(fun, x, N, sigma)

%   fun - Given as either a string or an inline function, returns the minus
%         log target and its gradient
%   x   - Initial state
%   N   - Number of steps
%   M   - Number of steps for running adaptation

[mllkd,mlp] = feval(fun, x);

np = length(x);
out.mlp = zeros(1,N);
out.samples = zeros(np,N);

dt  = exp(sigma);
a   = 2*sqrt(2*dt)/(2+dt);  % compute beta
b   = (2-dt)/(2+dt);        % sqrt( 1-beta^2 )

acc_t = 0;
% start MCMC
for i = 1:N
    %
    r   = randn(np,1);
    xn  = b*x + a*r;
    [mllkdn,mlpn] = feval(fun, xn);
    alpha   = mllkd - mllkdn;
    
    if  log(rand) < alpha
        x = xn;
        mlp = mlpn;
        mllkd = mllkdn;
        acc_t = acc_t + 1;
    end
    
    if mod(i,100)==0, disp(i), end
    
    %
    out.samples(:,i) = x;
    out.mlp(i) = mlp+mllkd;
end

out.acc_rate = acc_t/N;

end

