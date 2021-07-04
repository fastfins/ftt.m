classdef LagrangepCDF < Lagrangep & onedCDF
    
    properties
        cheby Chebyshev2nd
        cdf_data
    end
    
    methods
        function obj = LagrangepCDF(poly, varargin)
            obj@Lagrangep(poly.order, poly.num_elems, poly.domain, 'ghost_size', poly.gs, 'bc', poly.bc);
            obj@onedCDF(varargin{:});
            
            % local CDF poly
            obj.cheby = lag2cheby(poly.order*2, poly.local.domain);
            % vectorise    
            obj.global2local = reshape(obj.global2local, [], 1);
            obj.nodes = obj.nodes(:);
            %
            obj.cdf_basis2node = eval_int_basis(obj.cheby, obj.cheby.nodes);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cheby = lag2cheby(order, domain)
            %
            %Define data structure that maps lagrange polynomials to
            %Chebyshev polynomials with preserved boundary values. This
            %function is mainly useful for the squared mode--in which the
            %order of the polynomial is doubled.
            %
            %%%%%
            %Inputs:
            %
            %domain:
            %  domain of the input polynomial
            %
            %n_nodes:
            %  number of nodes used in the transformation
            %
            %
            %%%%%
            %Output:
            %
            %def:
            %  A data structure contains:
            %
            %  nodes:       nodes used for the transformation from Lagrange to Chebyshev
            %               we need to evaluate the pdf function on these nodes
            %  num_nodes:   number of nodes used in the transformation
            %  domain:      domain of the transformation
            %  basis2node:  inverse of the vandermonde matrix, transform nodal values
            %               to coefficient of
            %  node2basis:  Vandermonde matrix
            %
            %Tiangang Cui, August, 2019
            
            n_nodes = order + 1;
            cheby = Chebyshev2nd(order, domain);
            
            if n_nodes < 3
                error('Must use more than three nodes')
            end
            
            tmp   = Chebyshev2nd(n_nodes-3, domain);
            cheby.nodes       = [cheby.domain(1); tmp.nodes(:); cheby.domain(2)];
            cheby.basis2node  = eval_basis(cheby.nodes);
            
            [L,U] = lu(vpa(cheby.basis2node));
            cheby.node2basis  = double(U\(L\eye(n_nodes)));
            cheby.num_nodes   = n_nodes;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function data = pdf2cdf(obj, pdf)
            if (sum(pdf(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pdf(:)<0))])
            end
            
            data.size = size(pdf,2);
            if data.size > 1
                data.pdf_left  = pdf(1,:);
                data.pdf_right = pdf(end,:);
                % 1st coord: local, 2nd coord: elems, 3rd: pdfs
                local_pdf   = reshape(pdf(obj.global2local,:), obj.local.num_nodes, obj.num_elems, data.size);
                % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
                local_pdf   = permute(local_pdf, [1,3,2]);
                
                data.coef   = reshape(obj.local.node2basis*reshape(local_pdf, obj.local.num_nodes, []), ...
                    obj.local.num_nodes, data.size, obj.num_elems);
                
                data.cdf_grid   = zeros(obj.num_elems+1, data.size);
                data.cdf_nodes  = zeros(obj.num_nodes, data.size);
                data.base       = zeros(obj.num_elems, data.size);
                switch obj.bc
                    case{'Dirichlet'}
                        data.cdf_grid(1,:) = 0.5*data.pdf_left*obj.gs^2;
                    otherwise
                        %bug line below
                        %data.cdf_grid(2,:) = data.pdf_left*obj.gs;
                        data.cdf_grid(1,:) = data.pdf_left*obj.gs;
                end
                ind = reshape(obj.global2local, [], obj.num_elems);
                for i = 1:obj.num_elems
                    tmp             = (obj.cdf_basis2node*obj.jac(i))*data.coef(:,:,i);
                    data.base(i,:)  = tmp(1,:);
                    data.cdf_nodes(ind(:,i),:)  = tmp - data.base(i,:) + data.cdf_grid(i,:);
                    data.cdf_grid(i+1,:)        = data.cdf_grid(i,:) + tmp(end,:) - data.base(i,:);
                end
                switch obj.bc
                    case{'Dirichlet'}
                        data.norm = data.cdf_grid(end,:) + 0.5*data.pdf_right*obj.gs^2;
                    otherwise
                        data.norm = data.cdf_grid(end,:) + data.pdf_right*obj.gs;
                end
                data.cdf_nodes = data.cdf_nodes./data.norm;
                data.cdf_grid  = data.cdf_grid ./data.norm;
            else
                data.pdf_left  = pdf(1);
                data.pdf_right = pdf(end);
                % 1st coord: local, 2nd coord: elems
                local_pdf   = reshape(pdf(obj.global2local), obj.local.num_nodes, obj.num_elems);
                data.coef   = obj.local.node2basis*local_pdf;
                %
                data.cdf_grid   = zeros(obj.num_elems+1, 1);
                data.cdf_nodes  = zeros(obj.num_nodes, 1);
                data.base       = zeros(obj.num_elems, 1);
                switch obj.bc
                    case{'Dirichlet'}
                        data.cdf_grid(1) = 0.5*data.pdf_left*obj.gs^2;
                    otherwise
                        data.cdf_grid(1) = data.pdf_left*obj.gs;
                end
                ind = reshape(obj.global2local, [], obj.num_elems);
                for i = 1:obj.num_elems
                    tmp             = (obj.cdf_basis2node*obj.jac(i))*data.coef(:,i);
                    data.base(i)    = tmp(1);
                    data.cdf_nodes(ind(:,i)) = tmp - data.base(i) + data.cdf_grid(i);
                    data.cdf_grid(i+1) = data.cdf_grid(i) + tmp(end) - data.base(i);
                end
                switch obj.bc
                    case{'Dirichlet'}
                        data.norm = data.cdf_grid(end) + 0.5*data.pdf_right*obj.gs^2;
                    otherwise
                        data.norm = data.cdf_grid(end) + data.pdf_right*obj.gs;
                end
                %data.cdf_nodes = data.cdf_nodes/data.norm;
                %data.cdf_grid  = data.cdf_grid /data.norm;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function F = eval_int_lagp_local(obj, data, ei, mask, x)
            domains = [reshape(obj.grid(ei),[],1), reshape(obj.grid(ei+1),[],1)];
            x2z = 0.5*(domains(:,2)-domains(:,1));
            mid = 0.5*(domains(:,2)+domains(:,1));
            x   = (x(:)-mid)./x2z;            
            if data.size == 1
                b   = eval_int_basis(obj.cheby, x); % rewrite
                tmp = reshape(sum(b.*data.coef(:,ei)',2), size(x));
                F   = tmp - reshape(data.base(ei), size(tmp)) + reshape(data.cdf_grid(ei), size(tmp));
                %F   = (tmp - reshape(data.base(ei), size(tmp)))./data.norm + reshape(data.cdf_grid(ei), size(tmp));
            else
                %ii  = reshape(1:data.size, size(ei));
                ii  = 1:data_size;
                pi  = ii(mask);
                j1  = pi + (ei-1)*data.size;
                coe = data.coef(:,j1);
                
                b   = eval_int_basis(obj.cheby, x); % rewrite
                tmp = reshape(sum(b.*coe',2), size(x));
                
                j2  = (pi-1)* obj.num_elems    + ei;
                j3  = (pi-1)*(obj.num_elems+1) + ei;
                F   = tmp - reshape(data.base(j2),size(tmp)) + reshape(data.cdf_grid(j3),size(tmp));
                %F   = (tmp-reshape(data.base(j2),size(tmp)))./reshape(data.norm(mask),size(tmp)) + reshape(data.cdf_grid(j3),size(tmp));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_int_lag(obj, data, r)
            if data.size > 1 && data.size ~= length(r)
                error('Error: dimenion mismatch')
            end
            [mxi, nxi] = size(r);
            r = reshape(r,[],1);
            z = zeros(size(r));
            % index to locate the element
            ei = sum(reshape(obj.grid,1,[]) < r,2)';
            mask1 = ei==0;
            if sum(mask1) > 0
                tmp = ( r(mask1) - obj.domain(1) )./obj.gs ;
                switch obj.bc
                    case{'Dirichlet'}
                        if data.size == 1
                            z(mask1) = 0.5*(tmp.^2)*data.pdf_left;
                        else
                            z(mask1) = 0.5*(tmp.^2).*data.pdf_left(mask1)';
                        end
                    otherwise
                        if data.size == 1
                            z(mask1) = tmp*data.pdf_left;
                        else
                            z(mask1) = tmp.*data.pdf_left(mask1)';
                        end
                end
            end
            mask2 = ei==(obj.num_elems+1);
            if sum(mask2) > 0
                tmp = 1 - ( r(mask2) - obj.grid(end) )./obj.gs;
                switch obj.bc
                    case{'Dirichlet'}
                        if data.size == 1
                            z(mask2) = 1 - 0.5*(tmp.^2)*data.pdf_right;
                        else
                            z(mask2) = 1 - 0.5*(tmp.^2).*data.pdf_right(mask2)';
                        end
                    otherwise
                        if data.size == 1
                            z(mask2) = 1 - tmp*data.pdf_right;
                        else
                            z(mask2) = 1 - tmp.*data.pdf_right(mask2)';
                        end
                end
            end
            mask3 = ~mask1 & ~mask2;
            if sum(mask3) > 0
                z(mask3) = eval_int_lagp_local(obj, data, ei(mask3), mask3, r(mask3));
            end
            z = reshape(z, mxi, nxi);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf(obj, pk, r)
            data = pdf2cdf(obj, pk);
            z = eval_int_lag(obj, data, r);
            [mxi, nxi] = size(r);
            z = reshape(z(:)./data.norm(:), mxi, nxi);
        end
        
        %{
        function z = eval_cdf(obj, pk, r)
            data = pdf2cdf(obj, pk);
            if data.size > 1 && data.size ~= length(r)
                error('Error: dimenion mismatch')
            end
            [mxi, nxi] = size(r);
            r = reshape(r,[],1);
            z = zeros(size(r));
            % index to locate the element
            ei = sum(reshape(obj.grid,1,[]) < r,2)';
            mask1 = ei==0;
            if sum(mask1) > 0
                tmp = ( r(mask1) - obj.domain(1) )./obj.gs ;
                switch obj.bc
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
            mask2 = ei==(obj.num_elems+1);
            if sum(mask2) > 0
                tmp = 1 - ( r(mask2) - obj.grid(end) )./obj.gs;
                switch obj.bc
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
                if data.size == 1
                else
                    ii  = 1:data_size;
                    j1  = ii(mask3) + (ei(mask3)-1)*data.size;
                    coe = data.coef(:,j1);
                    z(mask3) = eval_cdf_local(obj, data, ei, mask3, r(mask3));
                end
            end
            z = reshape(z, mxi, nxi);
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf_deri(obj, pk, r)
            data = pdf2cdf(obj, pk);
            z = eval_int_lag(obj, data, r);
            [mxi, nxi] = size(r);
            z = reshape(z, mxi, nxi);
        end
        
        %{
        function z = eval_cdf_deri(obj, pk, r)
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
            
            if data_size == 1
                pdf_left  = pk(1);
                pdf_right = pk(end);
                % 1st coord: local, 2nd coord: elems
                local_pdf = reshape(pk(obj.global2local), obj.local.num_nodes, obj.num_elems);
                coef = obj.local.node2basis*local_pdf;
                %
                base = zeros(obj.num_elems, 1);
                cdf_grid = zeros(obj.num_elems+1, 1);
                switch obj.bc
                    case{'Dirichlet'}
                        cdf_grid(1) = 0.5*pdf_left*obj.gs^2;
                    otherwise
                        cdf_grid(1) = pdf_left*obj.gs;
                end
                for i = 1:obj.num_elems
                    base(i)  = (obj.cdf_basis2node(1,:)*obj.jac(i))*coef(:,i);
                    tmp_right = (obj.cdf_basis2node(end,:)*obj.jac(i))*coef(:,i);
                    cdf_grid(i+1) = cdf_grid(i) + tmp_right - base(i);
                end
            else
                pdf_left  = pk(1,:);
                pdf_right = pk(end,:);
                % 1st coord: local, 2nd coord: elems, 3rd: pdfs
                local_pdf = reshape(pk(obj.global2local,:), obj.local.num_nodes, obj.num_elems, []);
                % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
                local_pdf = permute(local_pdf, [1,3,2]);
                %
                coef = reshape(obj.local.node2basis*reshape(local_pdf, obj.local.num_nodes, []), ...
                    obj.local.num_nodes, [], obj.num_elems);
                %
                base = zeros(obj.num_elems, data_size);
                cdf_grid  = zeros(obj.num_elems+1, data_size);
                switch obj.bc
                    case{'Dirichlet'}
                        cdf_grid(1,:) = 0.5*pdf_left*obj.gs^2;
                    otherwise
                        cdf_grid(1,:) = pdf_left*obj.gs;
                end
                for i = 1:obj.num_elems
                    base(i,:) = (obj.cdf_basis2node(1,:)*obj.jac(i))*coef(:,:,i);
                    tmp_right = (obj.cdf_basis2node(end,:)*obj.jac(i))*coef(:,:,i);
                    cdf_grid(i+1,:) = cdf_grid(i,:) + tmp_right - base(i,:);
                end
            end
        
            %%%%%%%%%%%%%%%%%%%%
            
            z = zeros(size(r));
            % index to locate the element
            ei = sum(reshape(obj.grid,1,[]) < r,2)';
            
            mask1 = ei==0; % left ghost cell
            if sum(mask1) > 0
                tmp = ( r(mask1) - obj.domain(1) )./obj.gs ;
                switch obj.bc
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
            mask2 = ei==(obj.num_elems+1); % right ghost cell
            if sum(mask2) > 0
                tmp = 1 - ( r(mask2) - obj.grid(end) )./obj.gs;
                switch obj.bc
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
                %z = eval_lagrange_cdf_local(obj, data, ei, mask3, r(mask3));
                domains = [reshape(obj.grid(ei(mask3)),[],1), reshape(obj.grid(1+ei(mask3)),[],1)];
                if data_size == 1
                    %
                    b   = eval_oned_int_basis(obj.cheby, domains, obj.order, r(mask3));
                    tmp = reshape(sum(b.*coef(:,ei(mask3))',2), [], 1);
                    z   = (tmp - reshape(base(ei(mask3)), size(tmp))) + reshape(cdf_grid(ei(mask3)), size(tmp));
                else
                    %ii  = reshape(1:data_size, size(ei));
                    ii  = 1:data_size;
                    j1  = ii(mask3) + (ei(mask3)-1)*data_size;
                    coef = reshape(coef, obj.local.num_nodes, []);
                    %
                    b   = eval_oned_int_basis(obj.cheby, domains, obj.order, r(mask3));
                    tmp = reshape(sum(b.*coef(:,j1)',2), [], 1);
                    %
                    j2  = (ii(mask3)-1)* obj.num_elems    + ei(mask3);
                    j3  = (ii(mask3)-1)*(obj.num_elems+1) + ei(mask3);
                    z(mask3) = (tmp-reshape(base(j2),size(tmp))) + reshape(cdf_grid(j3),size(tmp));
                end
            end
            %{

            %}
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf(obj, pk, xi)
            
            data = pdf2cdf(obj, pk);
            
            if data.size > 1 && data.size ~= length(xi)
                error('Error: dimenion mismatch')
            end
            
            [mxi, nxi] = size(xi);
            xi = reshape(xi,[],1);
            
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
                switch obj.bc
                    case{'Dirichlet'}
                        r(mask1) = sqrt( 2*tmp )*obj.gs + obj.domain(1);
                    otherwise
                        r(mask1) = tmp*obj.gs + obj.domain(1);
                end
            end
            mask2 = ei==(obj.num_elems+1);
            if sum(mask2) > 0
                if data.size == 1
                    tmp = (1 - xi(mask2))*data.norm/data.pdf_right;
                else
                    tmp = (1 - xi(mask2)).*data.norm(mask2)'./data.pdf_right(mask2)';
                end
                tmp(isinf(tmp)) = 0;
                tmp(isnan(tmp)) = 0;
                tmp(tmp < eps)  = 0;
                switch obj.bc
                    case{'Dirichlet'}
                        r(mask2) = obj.domain(2) - sqrt(2*tmp)*obj.gs;
                    otherwise
                        r(mask2) = obj.domain(2) - tmp*obj.gs;
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
                a = obj.nodes(ind);
                b = obj.nodes(ind+1);
                %
                if data.size == 1
                rhs = data.norm*reshape(xi(mask3),[],1);
                else
                rhs = reshape(data.norm(mask3),[],1) .* reshape(xi(mask3),[],1);
                end
                r(mask3) = regula_falsi(@(x) eval_int_lagp_local(obj, data, ei(mask3), mask3, x) - rhs, obj.tol, a(:), b(:));
            end
            
            r = reshape(r, mxi, nxi);
            
            r(isnan(r)) = 0.5*(obj.domain(1) + obj.domain(2));
            r(isinf(r)) = 0.5*(obj.domain(1) + obj.domain(2));
        end
        
    end
    
end
