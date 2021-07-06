classdef spectralCDF < onedCDF
    
    properties
        sampling_nodes(1,:)
        cdf_basis2node(:,:)
    end
    
    methods (Abstract)
        eval_ref_int_basis(obj)
    end
        
    methods
        function obj = spectralCDF(varargin)
            obj@onedCDF(varargin{:});
            %
            obj.sampling_nodes      = linspace(obj.domain(1), obj.domain(2), max(obj.num_nodes*2, 200));
            obj.sampling_nodes(1)   = obj.domain(1)-eps;
            obj.sampling_nodes(end) = obj.domain(2)+eps;
            [x,J] = domain2reference(obj,obj.sampling_nodes(:));
            b = eval_ref_int_basis(obj, x);
            obj.cdf_basis2node      = b.*J;
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_int(obj, coef, x)
            [x,J] = domain2reference(obj,x(:));
            b = eval_ref_int_basis(obj, x);
            b = b.*J;
            %
            if size(coef,2) > 1
                if size(coef,2) == length(x)
                    f = reshape(sum(b.*coef',2), [], 1);
                else
                    error('Error: dimenion mismatch')
                end
            else
                f = reshape((b*coef), [], 1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf(obj, pk, r)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            if size(pk,2) > 1 && size(pk,2) ~= length(r)
                error('Error: dimenion mismatch')
            end
            r = reshape(r,[],1);
            %
            coef = obj.node2basis*pk;
            base = obj.cdf_basis2node(1,:)*coef;
            norm = (obj.cdf_basis2node(end,:)-obj.cdf_basis2node(1,:))*coef;
            %
            if size(pk,2) == 1
                z = eval_int(obj, coef, r) - base;
            else
                tmp = eval_int(obj, coef, r);
                z = reshape(tmp, size(r)) - reshape(base, size(r));
            end
            z = reshape(z(:)./norm(:), size(r));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf_deri(obj, pk, r)
            r = reshape(r,[],1);
            coef = obj.node2basis*pk;
            base = obj.cdf_basis2node(1,:)*coef;
            if size(pk,2) == 1
                z = eval_int(obj, coef, r) - base;
                z = reshape(z, size(r));
            else
                tmp = eval_int(obj, coef, r);
                z = (reshape(tmp, size(r)) - reshape(base, size(r)));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf(obj, pk, xi)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            data_size = size(pk,2);
            coef = obj.node2basis*pk;
            cdf_nodes = obj.cdf_basis2node*coef;
            cdf_base  = cdf_nodes(1,:);
            cdf_nodes = cdf_nodes - cdf_base;
            cdf_norm  = cdf_nodes(end,:);
            if data_size > 1 && data_size ~= length(xi)
                error('Error: dimenion mismatch')
            end
            %
            xi = reshape(xi,[],1);
            r = zeros(size(xi));
            
            if data_size == 1
                rhs = xi.*cdf_norm; % vertical
                ind = sum(reshape(cdf_nodes,1,[]) < rhs(:),2)';
            else
                rhs = xi(:).*cdf_norm(:); % vertical
                ind = sum(cdf_nodes < reshape(rhs,1,[]), 1);
            end
            mask1 = ind==0 | reshape(xi,1,[])<=eps;
            mask3 = ind==length(obj.sampling_nodes) | reshape(xi,1,[])>=1-eps;
            mask2 = ~(mask1|mask3);
            %
            r(mask1) = obj.sampling_nodes(1);
            r(mask3) = obj.sampling_nodes(end);
            %
            a = obj.sampling_nodes(ind(mask2));
            b = obj.sampling_nodes(ind(mask2)+1);
            %
            if data_size == 1
                r(mask2) = regula_falsi(@(x) eval_int(obj,coef,x)-(cdf_base+rhs(mask2)),...
                    obj.tol, a(:), b(:));
            else
                r(mask2) = regula_falsi(@(x) eval_int(obj,coef(:,mask2),x)-(reshape(cdf_base(mask2),[],1)+rhs(mask2)),...
                    obj.tol, a(:), b(:));
            end
            %
            r = reshape(r, size(xi));
            r(isnan(r)) = 0.5*(obj.domain(1) + obj.domain(2));
            r(isinf(r)) = 0.5*(obj.domain(1) + obj.domain(2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        function F = invert_obj(obj, data, mask, x)
            if data.size == 1
                F = eval_int(obj, coef(:,mask), x) - cdf_base(mask);
            else
                tmp = eval_int(obj, coef(:,mask), x);
                F = reshape(tmp, size(x)) - reshape(cdf_base(mask), size(x));
            end    
        end
        %}
    end
    
end