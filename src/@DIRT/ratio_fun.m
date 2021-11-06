function f = ratio_fun(obj, func, z, sqrt_flag)
% Evaluate the ratio function during the DIRT construction
%   f = RATIO_FUN(dirt, func, k, z, flag)
%
%   func - function handle of the target function, returns minus log
%          likelihood and minus log prior
%   flag - if return the sqrt of the target function, if FTT.sqrt_flag is
%          true, then this should be false
%   z    - reference random variables, d x n
%   f    - function value

[x, logft] = eval_irt(obj, z);

% compute the minus log likelihood and minus log prior
[mllkds, mlps] = func(x);

% compute the reference density at z
logfz = log_pdf(obj.diag, z);

if obj.n_layers > 0 % not in the first layer
    beta_p = obj.betas(obj.n_layers);
    beta = obj.betas(obj.n_layers+1);
    switch obj.method
        case {'Aratio'}
            logf =  (beta_p - beta)*mllkds + logfz;
        case {'Eratio'}
            logf =  - beta*mllkds - mlps - logft + logfz;
    end
else
    beta = obj.betas(obj.n_layers+1);
    logf = - beta*mllkds - mlps;
end

if sqrt_flag
    f = exp(0.5*logf);
else
    f = exp(logf);
end

end