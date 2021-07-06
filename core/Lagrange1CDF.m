classdef Lagrange1CDF < Lagrange1 & piecewiseCDF
    
    properties
        T
        iV
    end
    
    methods
        function obj = Lagrange1CDF(poly, varargin)
            obj@Lagrange1(poly.num_elems, poly.domain, 'ghost_size', poly.gs, 'bc', poly.bc);
            obj@piecewiseCDF(varargin{:});
            %
            obj.nodes = linspace(obj.domain(1)+obj.gs, obj.domain(2)-obj.gs, obj.num_elems*2+1);
            obj.num_nodes = length(obj.nodes);
            %
            dhh = obj.elem_size/2;
            %
            ii = zeros(3,obj.num_elems);
            jj = zeros(3,obj.num_elems);
            %
            % L = [1, 0, 0; 1, 1, 0; 1, 2, 1];
            % U = [1, 0, 0; 0, dhh, dhh^2; 0, 0, dhh^2*2];
            % LU = Vandermore
            % Vandermore = [1, 0, 0; 1, dhh, dhh^2; 1, dhh*2, 4*dhh^2];
            obj.iV = [1, 0, 0; -3/(2*dhh), 2/dhh, -1/(2*dhh); 1/(2*dhh^2), -1/dhh^2, 1/(2*dhh^2)];
            for i = 1:obj.num_elems
                ind = (1:3)+(i-1)*3;
                ii(:,i) = ind;
                jj(:,i) = (i-1)*2 + (1:3);
            end
            obj.T = sparse(ii(:), jj(:), ones(obj.num_elems*3,1), obj.num_elems*3, obj.num_nodes);
            %obj.TT = kron(speye(obj.num_elems),sparse(iV))*T;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function data = pdf2cdf(obj, pdf)
            data.size = size(pdf,2);
            if data.size > 1
                %{
                tmp = reshape(obj.TT*pdf, 3, []);
                data.A = reshape(tmp(1,:), obj.num_elems, []);
                data.B = reshape(tmp(2,:), obj.num_elems, []);
                data.C = reshape(tmp(3,:), obj.num_elems, []);
                %}
                
                % 1st coord: local, 2nd coord: elems, 3rd: pdfs
                data.coef = obj.iV*reshape(obj.T*pdf, 3, []);
                % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
                % data.coef = permute(reshape(obj.TT*pdf, 3, obj.num_elems, []), [1,3,2]);
                
                data.pdf_left  = pdf(1,:);
                data.pdf_right = pdf(end,:);
                [cdf_left, cdf_right] = pdf2cdf_bnd(obj, data.pdf_left, data.pdf_right);
                %
                %cdf_elems = data.A*obj.elem_size + data.B.*(obj.elem_size^2/2) + data.C.*(obj.elem_size^3/3);
                cdf_elems = reshape([obj.elem_size, obj.elem_size^2/2, obj.elem_size^3/3]*reshape(data.coef, 3, []), obj.num_elems, []);
                %
                data.cdf_grid  = zeros(obj.num_elems+1, data.size);
                data.cdf_grid(2:end,:) = cumsum(cdf_elems,1);
                data.cdf_grid = data.cdf_grid +  cdf_left;
                %
                data.norm = data.cdf_grid(end,:) + cdf_right;
                %
            else
                %{
                tmp = reshape(obj.TT*pdf, 3, []);
                data.A = reshape(tmp(1,:),[],1);
                data.B = reshape(tmp(2,:),[],1);
                data.C = reshape(tmp(3,:),[],1);
                %}
                data.coef = obj.iV*reshape(obj.T*pdf, 3, []);
                
                data.pdf_left  = pdf(1);
                data.pdf_right = pdf(end);
                [cdf_left, cdf_right] = pdf2cdf_bnd(obj, data.pdf_left, data.pdf_right);
                %
                %cdf_elems = data.A*obj.elem_size + data.B.*(obj.elem_size^2/2) + data.C.*(obj.elem_size^3/3);
                cdf_elems = [obj.elem_size, obj.elem_size^2/2, obj.elem_size^3/3]*reshape(data.coef, 3, []);
                %
                data.cdf_grid  = zeros(obj.num_elems+1, 1);
                data.cdf_grid(2:end) = cumsum(cdf_elems);
                data.cdf_grid = data.cdf_grid +  cdf_left;
                %
                data.norm = data.cdf_grid(end) + cdf_right;
                %
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_int_lag_local(obj, data, ei, mask, r)
            r = r - reshape(obj.grid(ei),size(r));
            if data.size > 1
                coi = find(mask);
                ind = ei + (coi-1)*obj.num_elems;
                jnd = ei + (coi-1)*(obj.num_elems+1);
                %
                % A: data.A(ind), B: data.B(ind), C: data.C(ind),
                % cdf: data.cdf_grid(jnd), base: data.base(ind)
                %
                % z = A(:).*r(:) + B(:).*r(:).^2/2 + C(:).*r(:).^3/3;
                % z = z + cdf(:);
                %
                f = sum([r(:), r(:).^2/2, r(:).^3/3].*data.coef(:,ind)', 2) + reshape(data.cdf_grid(jnd),[],1);
            else
                f = sum([r(:), r(:).^2/2, r(:).^3/3].*data.coef(:,ei)', 2) + reshape(data.cdf_grid(ei),[],1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [f,df] = eval_int_lag_local_newton(obj, data, ei, mask, rhs, r)
            r = r - reshape(obj.grid(ei),size(r));
            if data.size > 1
                coi = find(mask);
                ind = ei + (coi-1)*obj.num_elems;
                jnd = ei + (coi-1)*(obj.num_elems+1);
                %
                % A: data.A(ind), B: data.B(ind), C: data.C(ind),
                % cdf: data.cdf_grid(jnd), base: data.base(ind)
                %
                % z = A(:).*r(:) + B(:).*r(:).^2/2 + C(:).*r(:).^3/3;
                % z = z + cdf(:);
                %
                f  = sum([r(:), r(:).^2/2, r(:).^3/3].*data.coef(:,ind)', 2) + reshape(data.cdf_grid(jnd),[],1);
                df = sum([ones(length(r),1), r(:), r(:).^2].*data.coef(:,ind)', 2);
            else
                f  = sum([r(:), r(:).^2/2, r(:).^3/3].*data.coef(:,ei)', 2) + reshape(data.cdf_grid(ei),[],1);
                df = sum([ones(length(r),1), r(:), r(:).^2].*data.coef(:,ei)', 2);
            end
            f = f - rhs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf_local(obj, data, ei, mask, rhs)
            %
            a = obj.grid(ei);
            b = obj.grid(ei+1);
            %
            r = newton(obj, data, ei, mask, rhs(:), a(:), b(:));
            %if sum( r>b(:) | r<a(:) ) ~=0
            %    warning('newton failed')
            %    r = regula_falsi(obj, data, ei, mask, rhs(:), a(:), b(:));
            %end
        end

    end
    
end
