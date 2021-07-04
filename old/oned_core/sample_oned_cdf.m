function r = sample_oned_cdf(oned_cdf, pk, xi)
%Generating random variables using inverse transform. It supports two modes:
%  n cdfs and n random seeds
%  1 cdf  and n random seeds
%
%The root finding methods are modified from functions (regula_falsi and illinois)
%implemented in the chebfun system: https://www.chebfun.org/
%
%Tiangang Cui, August, 2019


data = pdf2cdf(oned_cdf, pk);

if data.size > 1 && data.size ~= length(xi)
    error('Error: dimenion mismatch')
end

[mxi, nxi] = size(xi);
xi = reshape(xi,[],1);

% single cdf, multiple random variable case
switch oned_cdf.type
    case{'Lagrange'}
        r = zeros(size(xi));
        % index to locate the element
        if data.size == 1
            ei = sum(reshape(data.cdf_grid,1,[]) < xi,2)';
        else
            ei = sum(data.cdf_grid <= reshape(xi,1,[]), 1);
        end
        mask1 = ei==0;
        if sum(mask1) > 0
            if data.size == 1
                tmp = xi(mask1)*data.norm/data.pdf_left;
            else
                tmp = xi(mask1).*data.norm(mask1)'./data.pdf_left(mask1)';
            end
            tmp(isinf(tmp)) = 0;
            tmp(isnan(tmp)) = 0;
            tmp(tmp < eps)  = 0;
            switch oned_cdf.bc
                case{'Dirichlet'}
                    r(mask1) = sqrt( 2*tmp )*oned_cdf.gs + oned_cdf.domain(1);
                otherwise
                    r(mask1) = tmp*oned_cdf.gs + oned_cdf.domain(1);
            end
        end
        mask2 = ei==(oned_cdf.num_elems+1);
        if sum(mask2) > 0                 
            if data.size == 1
                tmp = (1 - xi(mask2))*data.norm/data.pdf_right;
            else
                tmp = (1 - xi(mask2)).*data.norm(mask2)'./data.pdf_right(mask2)';
            end
            tmp(isinf(tmp)) = 0;
            tmp(isnan(tmp)) = 0;
            tmp(tmp < eps)  = 0;
            switch oned_cdf.bc
                case{'Dirichlet'}
                    %{
                    if data.size == 1
                        tmp = 1 - 2*(xi(mask2)-data.cdf_grid(end))*data.norm/data.pdf_right;
                    else
                        tmp = 1 - 2*(xi(mask2)-data.cdf_grid(end,mask2)').*data.norm(mask2)'./data.pdf_right(mask2)';
                    end
                    tmp(isinf(tmp)) = 0;
                    tmp(isnan(tmp)) = 0;
                    tmp(tmp < eps)  = 0;
                    r(mask2) = (1 - sqrt(tmp))*oned_cdf.gs + oned_cdf.elem_right(oned_cdf.num_elems);
                    %}

                    r(mask2) = oned_cdf.domain(2) - sqrt(2*tmp)*oned_cdf.gs;
                otherwise
                    r(mask2) = oned_cdf.domain(2) - tmp*oned_cdf.gs;
            end
        end
        
        mask3 = ~mask1 & ~mask2;
        if sum(mask3) > 0
            % index to locate the roots, the roots are located between ind
            % and ind+1 nodes
            if data.size == 1
                ind = sum(reshape(data.cdf_nodes,1,[]) <= reshape(xi(mask3),[],1), 2)';
            else
                ind = sum(data.cdf_nodes(:,mask3) <= reshape(xi(mask3),1,[]), 1);
            end
            a = oned_cdf.nodes(ind);
            b = oned_cdf.nodes(ind+1);
            r(mask3) = root_finding(@(x) eval_lagrange_cdf_local(oned_cdf, data, ei, mask3, x) - xi(mask3), ...
                oned_cdf.err_tol, oned_cdf.method, a(:), b(:));
        end
        
    otherwise
        % 'Chebyshev1st', 'Chebyshev2nd','Legendre', 'Jacobi11', 'Fourier'
        % index to locate the roots, the roots are located between ind
        % and ind+1 nodes
        if data.size == 1
            ind = sum(reshape(data.cdf_nodes,1,[]) <= xi(:),2)';
        else
            ind = sum(data.cdf_nodes <= reshape(xi,1,[]), 1);
        end
        ind(ind==0) = 1;
        a = oned_cdf.sampling_nodes(ind);
        b = oned_cdf.sampling_nodes(ind+1);
        r = root_finding(@(x) eval_orthogonal_cdf(oned_cdf, data, x) - xi, oned_cdf.err_tol, oned_cdf.method, a(:), b(:));
end

r = reshape(r, mxi, nxi);

r(isnan(r)) = 0.5*(oned_cdf.domain(1) + oned_cdf.domain(2));
r(isinf(r)) = 0.5*(oned_cdf.domain(1) + oned_cdf.domain(2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = root_finding(func, tol, method, a, b)

switch method
    case{'Regula_Falsi'}
        r = regula_falsi(func, tol, a, b);
    case{'Illinois'}
        r = illinois(func, tol, a, b);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = regula_falsi(func, tol, a, b)


fa = func(a) ;
fb = func(b) ;

if sum(sign(fb.*fa) ~= -1)
    disp('Root finding: initial guesses on one side')
end

c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cold    = inf;

while ( norm(c-cold, inf) >= tol )
    cold    = c;
    fc  = func(c) ;
    
    I1  = (fc < 0);
    I2  = (fc > 0);
    I3  = ~I1 & ~I2;
    a   = I1.*c + I2.*a + I3.*c;
    b   = I1.*b + I2.*c + I3.*c;
    fa  = I1.*fc + I2.*fa + I3.*fc;
    fb  = I1.*fb + I2.*fc + I3.*fc;
    step    = -fb.*(b - a)./(fb - fa);
    step(isnan(step)) = 0;
    c = b + step;
    
    %norm(fc, inf)
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = illinois(func, tol, a, b)

fa = func(a) ;
fb = func(b) ;

c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cold    = inf;

side = zeros(size(a));

while ( norm(c - cold, inf) >= tol )
    cold    = c;
    fc  = func(c) ;
    
    I1  = (fc < 0);
    I2  = (fc > 0);
    I3  = ~I1 & ~I2;
    a   = I1.*c + I2.*a + I3.*c;
    b   = I1.*b + I2.*c + I3.*c;
    fa  = I1.*fc + I2.*fa + I3.*fc;
    fb  = I1.*fb + I2.*fc + I3.*fc;
    
    fb(side(I1) == -1) = fb(side(I1) == -1)/2;
    side(I1) = -1;
    fa(side(I1) == 1) = fa(side(I1) == 1)/2;
    side(I2) = 1;
    
    step = -fb.*(b - a)./(fb - fa);
    step(isnan(step)) = 0;
    c = b + step;
end

end

