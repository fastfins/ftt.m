function f = eval_ou_process_marginal(data, ind, x)

C = data.C(ind, ind);
f = exp(-0.5*sum((C\x).*x,1));
z = sqrt(det(C)*(2*pi)^length(ind));
f = f/z;

end