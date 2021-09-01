function f = ratio_fun(obj, func, k, z, sqrt_flag)
% Evaluate the ratio function during the DIRT construction
%   f = RATIO_FUN(dirt, func, tp, t, z)
%
%   func - function handle of the target function, returns minus log
%          likelihood and minus log prior
%   tp   - previous temperature
%   t    - current temperature
%   z    - reference random variables, d x n
%   f    - function value

[x, logft] = eval_irt(obj, z);

% compute the minus log likelihood and minus log prior
[mllkds, mlps] = func(x);

% compute the reference density at z
logfz = log_pdf(obj.diag, z);

if k > 1
    beta_p = obj.betas(k-1);
    beta = obj.betas(k);
    switch obj.method
        case {'Aratio'}
            logf =  (beta_p - beta)*mllkds + logfz;
        case {'Eratio'}
            logf =  - beta*mllkds - mlps - logft + logfz;
    end
else
    beta = obj.betas(k);
    logf = - beta*mllkds - mlps;
end

if sqrt_flag
    f = exp(0.5*logf);
else
    f = exp(logf);
end

end