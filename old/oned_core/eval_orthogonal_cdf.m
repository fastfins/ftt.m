function F = eval_orthogonal_cdf(oned_cdf, data, x)
%Evaluate one dimensional cdf defined by orthogonal polynomials
%
%Tiangang Cui, August, 2019

if data.size == 1
    F = eval_oned_int(oned_cdf, data.coef, x) - data.base;
    F = F/data.norm;
else
    tmp = eval_oned_int(oned_cdf, data.coef, x);
    F   = (reshape(tmp, size(x)) - reshape(data.base, size(x)))./reshape(data.norm, size(x));
end

end