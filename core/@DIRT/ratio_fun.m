function f = ratio_fun(obj, func, beta_p, beta, z)
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
logfz = eval_log_pdf(obj.diag, z);
            
switch obj.method
    case {'Aratio'}
        f =  exp( (beta_p - beta)*mllkds - mlps + logfz);
    case {'Eratio'}
        f =  exp( - beta*mllkds - mlps - logft + logfz);
end

end