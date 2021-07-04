function x = sample_oned_domain(oned, n)
%
%Draw n uniform random variables from the domain of the polynomial defined 
%by the input oned. 
%
%Output: d x n matrix, where n is the number of random variables
%
%Tiangang Cui, August, 2019

switch oned.type
    case{'Lagrange'}
        x = rand(1,n)*(oned.elem_right(end)-oned.elem_left(1))+oned.elem_left(1);
    case{'Chebyshev1st','Chebyshev2nd', 'Legendre', 'Jacobi11', 'Fourier'}
        x = rand(1,n)*(oned.domain(2)-oned.domain(1))+oned.domain(1);
    case{'Laguerre'}
    case{'Hermite'}
end



end