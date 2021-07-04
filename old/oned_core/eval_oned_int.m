function f = eval_oned_int(oned_cdf, coef, x)
%Evaluate the integral of polynomial defined by chebyshev basis functions 
%at input x.
%
%Tiangang Cui, August, 2019

b = eval_oned_int_basis(oned_cdf.type, oned_cdf.domain, oned_cdf.order, x);

if size(coef,2) > 1
    if size(coef,2) == length(x)
        f = reshape(sum(b.*coef',2), size(x));
    else
        error('Error: dimenion mismatch')
    end
else
    f = reshape((b*coef), size(x));
end

end
% use the following normalisation
%
%
%
%normalising2nd = sqrt(2/pi);