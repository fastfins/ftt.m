function z = eval_oned_cdf_deri(oned_cdf, pk, r)
%Generating random variables using inverse transform. It supports two modes:
%  n cdfs and n random seeds
%  1 cdf  and n random seeds
%
%The root finding methods are modified from functions (regula_falsi and illinois)
%implemented in the chebfun system: https://www.chebfun.org/
%
%Tiangang Cui, August, 2019

r = reshape(r,[],1);
data_size = size(pk,2);

switch oned_cdf.type
    case{'Lagrange'}
        % compute coefficients
        if data_size == 1
            pdf_left  = pk(1);
            pdf_right = pk(end);
            % 1st coord: local, 2nd coord: elems
            local_pdf = reshape(pk(oned_cdf.global2local), oned_cdf.local.num_nodes, oned_cdf.num_elems);
            coef = oned_cdf.local.node2basis*local_pdf;
            %
            base = zeros(oned_cdf.num_elems, 1);
            cdf_grid = zeros(oned_cdf.num_elems+1, 1);
            switch oned_cdf.bc
                case{'Dirichlet'}
                    cdf_grid(1) = 0.5*pdf_left*oned_cdf.gs^2;
                otherwise
                    cdf_grid(1) = pdf_left*oned_cdf.gs;
            end
            for i = 1:oned_cdf.num_elems
                base(i)  = (oned_cdf.cdf_basis2node(1,:)*oned_cdf.jac(i))*coef(:,i);
                tmp_right = (oned_cdf.cdf_basis2node(end,:)*oned_cdf.jac(i))*coef(:,i);
                cdf_grid(i+1) = cdf_grid(i) + tmp_right - base(i);
            end
        else
            pdf_left  = pk(1,:);
            pdf_right = pk(end,:);
            % 1st coord: local, 2nd coord: elems, 3rd: pdfs
            local_pdf = reshape(pk(oned_cdf.global2local,:), oned_cdf.local.num_nodes, oned_cdf.num_elems, []);
            % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
            local_pdf = permute(local_pdf, [1,3,2]);
            %
            coef = reshape(oned_cdf.local.node2basis*reshape(local_pdf, oned_cdf.local.num_nodes, []), ...
                oned_cdf.local.num_nodes, [], oned_cdf.num_elems);
            %
            base = zeros(oned_cdf.num_elems, data_size);
            cdf_grid  = zeros(oned_cdf.num_elems+1, data_size);
            switch oned_cdf.bc
                case{'Dirichlet'}
                    cdf_grid(1,:) = 0.5*pdf_left*oned_cdf.gs^2;
                otherwise
                    cdf_grid(1,:) = pdf_left*oned_cdf.gs;
            end
            for i = 1:oned_cdf.num_elems
                base(i,:) = (oned_cdf.cdf_basis2node(1,:)*oned_cdf.jac(i))*coef(:,:,i);
                tmp_right = (oned_cdf.cdf_basis2node(end,:)*oned_cdf.jac(i))*coef(:,:,i);
                cdf_grid(i+1,:) = cdf_grid(i,:) + tmp_right - base(i,:);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%
        
        z = zeros(size(r));
        % index to locate the element
        ei = sum(reshape(oned_cdf.grid_pts,1,[]) < r,2)';
        
        mask1 = ei==0; % left ghost cell
        if sum(mask1) > 0
            tmp = ( r(mask1) - oned_cdf.domain(1) )./oned_cdf.gs ;
            switch oned_cdf.bc
                case{'Dirichlet'}
                    if data_size == 1
                        z(mask1) = 0.5*(tmp.^2)*pdf_left;
                    else
                        z(mask1) = 0.5*(tmp.^2).*pdf_left(mask1)';
                    end
                otherwise
                    if data_size == 1
                        z(mask1) = tmp*pdf_left;
                    else
                        z(mask1) = tmp.*pdf_left(mask1)';
                    end
            end
        end
        mask2 = ei==(oned_cdf.num_elems+1); % right ghost cell
        if sum(mask2) > 0
            tmp = 1 - ( r(mask2) - oned_cdf.elem_right(end) )./oned_cdf.gs;
            switch oned_cdf.bc
                case{'Dirichlet'}
                    if data_size == 1
                        z(mask2) = 1 - 0.5*(tmp.^2)*pdf_right;
                    else
                        z(mask2) = 1 - 0.5*(tmp.^2).*pdf_right(mask2)';
                    end
                otherwise
                    if data_size == 1
                        z(mask2) = 1 - tmp*pdf_right;
                    else
                        z(mask2) = 1 - tmp.*pdf_right(mask2)';
                    end
            end
        end
        
        mask3 = ~mask1 & ~mask2; % middle cell
        if sum(mask3) > 0
            %
            %z = eval_lagrange_cdf_local(oned_cdf, data, ei, mask3, r(mask3));
            domains = [reshape(oned_cdf.elem_left(ei(mask3)),[],1), reshape(oned_cdf.elem_right(ei(mask3)),[],1)];
            if data_size == 1
                %
                b   = eval_oned_int_basis(oned_cdf.cheby, domains, oned_cdf.order, r(mask3));
                tmp = reshape(sum(b.*coef(:,ei(mask3))',2), [], 1);
                z   = (tmp - reshape(base(ei(mask3)), size(tmp))) + reshape(cdf_grid(ei(mask3)), size(tmp));
            else
                ii  = reshape(1:data_size, size(ei));
                j1  = ii(mask3) + (ei(mask3)-1)*data_size;
                coef = reshape(coef, oned_cdf.local.num_nodes, []);
                %
                b   = eval_oned_int_basis(oned_cdf.cheby, domains, oned_cdf.order, r(mask3));
                tmp = reshape(sum(b.*coef(:,j1)',2), [], 1);
                %
                j2  = (ii(mask3)-1)* oned_cdf.num_elems    + ei(mask3);
                j3  = (ii(mask3)-1)*(oned_cdf.num_elems+1) + ei(mask3);
                z(mask3) = (tmp-reshape(base(j2),size(tmp))) + reshape(cdf_grid(j3),size(tmp));
            end
        end
        
    otherwise
        % 'Chebyshev1st', 'Chebyshev2nd','Legendre', 'Jacobi11', 'Fourier'
        coef = oned_cdf.node2basis*pk;
        cdf_nodes = oned_cdf.cdf_basis2node*coef;
        base = cdf_nodes(1,:);
        if data_size == 1
            z = eval_oned_int(oned_cdf, coef, r) - base;
        else
            tmp = eval_oned_int(oned_cdf, coef, r);
            z = (reshape(tmp, size(r)) - reshape(base, size(r)));
        end
end


