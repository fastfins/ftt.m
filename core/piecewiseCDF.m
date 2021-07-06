classdef piecewiseCDF < onedCDF
    
    methods (Abstract)
        eval_int_lag_local(obj)
        pdf2cdf(obj)
    end
    
    methods
        function obj = piecewiseCDF(varargin)
            obj@onedCDF(varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [cdf_left, cdf_right] = pdf2cdf_bnd(obj, pdf_left, pdf_right)
            switch obj.bc
                case{'Dirichlet'}
                    cdf_left  = 0.5*pdf_left*obj.gs^2;
                    cdf_right = 0.5*pdf_right*obj.gs^2;
                otherwise
                    cdf_left  = pdf_left*obj.gs;
                    cdf_right = pdf_right*obj.gs;
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
                            z(mask2) = data.norm - 0.5*(tmp.^2)*data.pdf_right;
                        else
                            z(mask2) = data.norm(mask2)' - 0.5*(tmp.^2).*data.pdf_right(mask2)';
                        end
                    otherwise
                        if data.size == 1
                            z(mask2) = data.norm - tmp*data.pdf_right;
                        else
                            z(mask2) = data.norm(mask2)' - tmp.*data.pdf_right(mask2)';
                        end
                end
            end
            mask3 = ~mask1 & ~mask2;
            if sum(mask3) > 0
                z(mask3) = eval_int_lag_local(obj, data, ei(mask3), mask3, r(mask3));
            end
            z = reshape(z, mxi, nxi);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf(obj, pk, r)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            data = pdf2cdf(obj, pk);
            z = eval_int_lag(obj, data, r);
            z = reshape(z(:)./data.norm(:), size(r));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf_deri(obj, pk, r)
            data = pdf2cdf(obj, pk);
            z = eval_int_lag(obj, data, r);
            z = reshape(z, size(r));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf(obj, pk, xi)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            data = pdf2cdf(obj, pk);
            if data.size > 1 && data.size ~= length(xi)
                error('Error: dimenion mismatch')
            end
            
            [mxi, nxi] = size(xi);
            xi = reshape(xi,[],1); % vertical
            
            r = zeros(size(xi)); % vertical
            % index to locate the element
            if data.size == 1
                rhs = xi.*data.norm;
                ei  = sum(reshape(data.cdf_grid,1,[]) < rhs,2)';
            else
                rhs = xi(:).*data.norm(:);
                ei  = sum(data.cdf_grid <= reshape(rhs,1,[]), 1);
            end
            mask1 = ei==0;
            if sum(mask1) > 0
                if data.size == 1
                    tmp = rhs(mask1)/data.pdf_left;
                else
                    tmp = rhs(mask1)./reshape(data.pdf_left(mask1),[],1);
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
                    tmp = (data.norm - rhs(mask2))/data.pdf_right;
                else
                    tmp = (reshape(data.norm(mask2),[],1) - rhs(mask2))./reshape(data.pdf_right(mask2),[],1);
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
            %
            mask3 = ~mask1 & ~mask2;
            if sum(mask3) > 0
                r(mask3) = invert_cdf_local(obj, data, ei(mask3), mask3, rhs(mask3));
            end
            %
            r = reshape(r, mxi, nxi);
            r(isnan(r)) = 0.5*(obj.domain(1) + obj.domain(2));
            r(isinf(r)) = 0.5*(obj.domain(1) + obj.domain(2));
        end
    end
end