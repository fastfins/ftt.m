function [u,f] = eval_irt(obj, v)


f = zeros(1, size(v,2));

u = v;
for l = obj.n_layers:-1:1
    z = 0.5*( 1 + erf(u/sqrt(2))); 
    %
    [u,fl] = eval_irt(obj.irts{l}, z);
    logref = - 0.5*sum(ur.^2, 1) - 0.5*log(2*pi)*size(ur,1);
    f = f + log(fr) - logref;
end

end