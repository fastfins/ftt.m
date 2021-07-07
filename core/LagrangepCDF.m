classdef LagrangepCDF < Lagrangep & piecewiseCDF
    
    properties
        cheby Chebyshev2nd
        cdf_basis2node(:,:)
    end
    
    methods (Static)
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
            cheby.basis2node  = eval_basis(cheby, cheby.nodes);
            
            [L,U] = lu(vpa(cheby.basis2node));
            cheby.node2basis  = double(U\(L\eye(n_nodes)));
            cheby.num_nodes   = n_nodes;
        end
    end
    
    
    methods
        function obj = LagrangepCDF(poly, varargin)
            obj@Lagrangep(poly.order, poly.num_elems, poly.domain, 'ghost_size', poly.gs, 'bc', poly.bc);
            obj@piecewiseCDF(varargin{:});
            
            % local CDF poly
            obj.cheby = LagrangepCDF.lag2cheby(poly.order*2, poly.local.domain);
            %
            obj.mass    = [];
            obj.mass_R  = [];
            obj.int_W   = [];
            obj.weights = [];
            %
            % setup global nodes
            obj.num_nodes   = obj.num_elems*(obj.cheby.num_nodes-1)+1;
            obj.nodes       = zeros(1, obj.num_nodes);
            for i = 1:obj.num_elems
                ind = ( 1:obj.cheby.num_nodes ) + (obj.cheby.num_nodes-1)*(i-1);
                obj.nodes(ind) = obj.cheby.nodes*obj.elem_size + obj.grid(i);
            end
            
            % map the function value y to each local element
            if obj.cheby.num_nodes > 2
                j = obj.cheby.num_nodes:(obj.cheby.num_nodes-1):obj.num_nodes;
                obj.global2local = reshape([reshape(1:(obj.num_nodes-1), obj.cheby.num_nodes-1, obj.num_elems); j], [], 1);
            else
                obj.global2local = reshape([1:(obj.num_nodes-1); 2:obj.num_nodes], [], 1);
            end
            % vectorise
            obj.global2local = reshape(obj.global2local, [], 1);
            obj.nodes = obj.nodes(:);
            %
            [x,J] = domain2reference(obj.cheby,obj.cheby.nodes(:));
            b = eval_ref_int_basis(obj.cheby, x);
            obj.cdf_basis2node = b.*J;
        end       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function data = pdf2cdf(obj, pdf)
            data.size = size(pdf,2);
            if data.size > 1
                data.pdf_left  = pdf(1,:);
                data.pdf_right = pdf(end,:);
                % 1st coord: local, 2nd coord: elems, 3rd: pdfs
                local_pdf   = reshape(pdf(obj.global2local,:), obj.cheby.num_nodes, obj.num_elems, data.size);
                % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
                local_pdf   = permute(local_pdf, [1,3,2]);
                
                data.coef   = reshape(obj.cheby.node2basis*reshape(local_pdf, obj.cheby.num_nodes, []), ...
                    obj.cheby.num_nodes, data.size, obj.num_elems);
                
                [cdf_left, cdf_right] = pdf2cdf_bnd(obj, data.pdf_left, data.pdf_right);
                data.cdf_grid   = zeros(obj.num_elems+1, data.size);
                data.cdf_nodes  = zeros(obj.num_nodes, data.size);
                data.base       = zeros(obj.num_elems, data.size);
                %
                data.cdf_grid(1,:) = cdf_left;
                %
                ind = reshape(obj.global2local, [], obj.num_elems);
                for i = 1:obj.num_elems
                    tmp             = (obj.cdf_basis2node*data.coef(:,:,i))*obj.jac;
                    data.base(i,:)  = tmp(1,:);
                    data.cdf_nodes(ind(:,i),:)  = tmp - data.base(i,:) + data.cdf_grid(i,:);
                    data.cdf_grid(i+1,:)        = data.cdf_grid(i,:) + tmp(end,:) - data.base(i,:);
                end
                %
                data.norm = data.cdf_grid(end,:) + cdf_right;
            else
                data.pdf_left  = pdf(1);
                data.pdf_right = pdf(end);
                % 1st coord: local, 2nd coord: elems
                local_pdf   = reshape(pdf(obj.global2local), obj.cheby.num_nodes, obj.num_elems);
                data.coef   = obj.cheby.node2basis*local_pdf;
                %
                [cdf_left, cdf_right] = pdf2cdf_bnd(obj, data.pdf_left, data.pdf_right);
                %
                data.cdf_grid   = zeros(obj.num_elems+1, 1);
                data.cdf_nodes  = zeros(obj.num_nodes, 1);
                data.base       = zeros(obj.num_elems, 1);
                %
                data.cdf_grid(1,:) = cdf_left;
                %
                ind = reshape(obj.global2local, [], obj.num_elems);
                for i = 1:obj.num_elems
                    tmp             = (obj.cdf_basis2node*data.coef(:,i))*obj.jac;
                    data.base(i)    = tmp(1);
                    data.cdf_nodes(ind(:,i)) = tmp - data.base(i) + data.cdf_grid(i);
                    data.cdf_grid(i+1) = data.cdf_grid(i) + tmp(end) - data.base(i);
                end
                %
                data.norm = data.cdf_grid(end,:) + cdf_right;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function F = eval_int_lag_local(obj, data, ei, mask, x)
            domains = [reshape(obj.grid(ei),[],1), reshape(obj.grid(ei+1),[],1)];
            x2z = 0.5*(domains(:,2)-domains(:,1));
            mid = 0.5*(domains(:,2)+domains(:,1));
            x   = (x(:)-mid)./x2z;
            if data.size == 1
                b   = eval_ref_int_basis(obj.cheby, x); % rewrite
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.coef(:,ei)',2), size(x));
                F   = tmp - reshape(data.base(ei),size(tmp)) + reshape(data.cdf_grid(ei),size(tmp));
            else
                pi  = find(mask);
                j1  = pi + (ei-1)*data.size;
                %
                b   = eval_ref_int_basis(obj.cheby, x); % rewrite
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.coef(:,j1)',2), size(x));
                %
                j2  = (pi-1)* obj.num_elems    + ei;
                j3  = (pi-1)*(obj.num_elems+1) + ei;
                F   = tmp - reshape(data.base(j2),size(tmp)) + reshape(data.cdf_grid(j3),size(tmp));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [f,df] = eval_int_lag_local_newton(obj, data, ei, mask, rhs, x)
            domains = [reshape(obj.grid(ei),[],1), reshape(obj.grid(ei+1),[],1)];
            x2z = 0.5*(domains(:,2)-domains(:,1));
            mid = 0.5*(domains(:,2)+domains(:,1));
            x   = (x(:)-mid)./x2z;
            if data.size == 1
                %b  = eval_ref_int_basis(obj.cheby, x); % rewrite
                [b,db] = eval_ref_int_basis_newton(obj.cheby, x);
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.coef(:,ei)',2), size(x));
                f   = tmp - reshape(data.base(ei),size(tmp)) + reshape(data.cdf_grid(ei),size(tmp));
                df  = reshape(sum(db.*data.coef(:,ei)',2), size(x));
            else
                pi  = find(mask);
                j1  = pi + (ei-1)*data.size;
                %
                %b  = eval_ref_int_basis(obj.cheby, x); % rewrite
                [b,db] = eval_ref_int_basis_newton(obj.cheby, x);
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.coef(:,j1)',2), size(x));
                %
                j2  = (pi-1)* obj.num_elems    + ei;
                j3  = (pi-1)*(obj.num_elems+1) + ei;
                f   = tmp - reshape(data.base(j2),size(tmp)) + reshape(data.cdf_grid(j3),size(tmp));
                df  = reshape(sum(db.*data.coef(:,j1)',2), size(x));
            end
            f = f - rhs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf_local(obj, data, ei, mask, rhs)
            if data.size == 1
                ind = sum(reshape(data.cdf_nodes,1,[]) <= reshape(rhs,[],1), 2)';
            else
                ind = sum(data.cdf_nodes(:,mask) <= reshape(rhs,1,[]), 1);
            end
            a = obj.nodes(ind);
            b = obj.nodes(ind+1);
            %
            %r = regula_falsi(obj, data, ei, mask, rhs(:), a(:), b(:));
            r = newton(obj, data, ei, mask, rhs(:), a(:), b(:));
        end
        
    end
    
end
