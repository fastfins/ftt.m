function z = eval_oned_cdf(oned_cdf, pk, r)
%Generating random variables using inverse transform. It supports two modes:
%  n cdfs and n random seeds
%  1 cdf  and n random seeds
%
%The root finding methods are modified from functions (regula_falsi and illinois)
%implemented in the chebfun system: https://www.chebfun.org/
%
%Tiangang Cui, August, 2019


data = pdf2cdf(oned_cdf, pk);

if data.size > 1 && data.size ~= length(r)
    error('Error: dimenion mismatch')
end

[mxi, nxi] = size(r);
r = reshape(r,[],1);

% single cdf, multiple random variable case
switch oned_cdf.type
    case{'Lagrange'}
        z = zeros(size(r));
        % index to locate the element
        ei = sum(reshape(data.grid_pts,1,[]) < r,2)';
                
        mask1 = ei==0;
        if sum(mask1) > 0
            tmp = ( r(mask1) - oned_cdf.domain(1) )./oned_cdf.gs ;
            switch oned_cdf.bc
                case{'Dirichlet'}
                    if data.size == 1
                        z(mask1) = 0.5*(tmp.^2)*data.pdf_left/data.norm;
                    else
                        z(mask1) = 0.5*(tmp.^2).*data.pdf_left(mask1)'./data.norm(mask1)';
                    end
                otherwise
                    if data.size == 1
                        z(mask1) = tmp*data.pdf_left/data.norm;
                    else
                        z(mask1) = tmp.*data.pdf_left(mask1)'./data.norm(mask1)';
                    end
            end
        end
        mask2 = ei==(oned_cdf.num_elems+1);
        if sum(mask2) > 0
            tmp = 1 - ( r(mask2) - oned_cdf.elem_right(end) )./oned_cdf.gs;
            switch oned_cdf.bc
                case{'Dirichlet'}
                    if data.size == 1
                        z(mask2) = 1 - 0.5*(tmp.^2)*data.pdf_right/data.norm;
                    else
                        z(mask2) = 1 - 0.5*(tmp.^2).*data.pdf_right(mask2)'./data.norm(mask2)';
                    end
                otherwise
                    if data.size == 1
                        z(mask2) = 1 - tmp*data.pdf_right/data.norm;
                    else
                        z(mask2) = 1 - tmp.*data.pdf_right(mask2)'./data.norm(mask2)';
                    end
            end
        end
        
        mask3 = ~mask1 & ~mask2;
        if sum(mask3) > 0
            z(mask3) = eval_lagrange_cdf_local(oned_cdf, data, ei, mask3, r(mask3));
        end
        
    otherwise
        % 'Chebyshev1st', 'Chebyshev2nd','Legendre', 'Jacobi11', 'Fourier'
        z = eval_orthogonal_cdf(oned_cdf, data, r);
end

z = reshape(z, mxi, nxi);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

