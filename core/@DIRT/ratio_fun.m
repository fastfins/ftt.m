function f = ratio_fun(obj, func, beta_p, beta, v)

[u, f_dirt] = eval_dirt(obj, v);

% compute the minus log likelihood and minus log prio
[mllkds, mlps] = func(u);

switch obj.method
    case {'Aratio'}
        f =  exp( (beta_p - beta)*mllkds - mlps );
    case {'Eratio'}
        f =  exp( - beta*mllkds - mlps ) / f_dirt;
end

end