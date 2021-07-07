classdef piecewiseCDF < onedCDF
    % piecewiseCDF class 
    %
    % For piecewise linear basis (Lagrange1), we directly apply Newton's
    % method after a grid search based on the Lagrange nodes. 
    %
    % For piecewise high order basis (Lagrangep), we first convert each 
    % piesewise Lagerange polynomial into 2nd Chebyshev basis, and then
    % work out the CDF. Before applying root findings, a grid search based 
    % on Chebyshev nodes is applied to locate the left and right boundary 
    % of root finding. The Chebyshev nodes contains both ends of the
    % interpolation interval. Newton's method is fragile towards the end of
    % the interval, so it's been modified to stay in the search interval. 
    %
    % See also Lagrange1CDF and LagrangepCDF.
    
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function c = regula_falsi(obj, data, ei, mask, rhs, a, b)    
            fa = eval_int_lag_local(obj, data, ei, mask, a) - rhs;
            fb = eval_int_lag_local(obj, data, ei, mask, b) - rhs;
            if sum(sign(fb.*fa) ~= -1)
                disp('Root finding: initial guesses on one side')
            end
            c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
            cold = inf;
            %i = 2;
            while ( norm(c-cold, Inf) > obj.tol )
                cold = c;
                fc  = eval_int_lag_local(obj, data, ei, mask, c) - rhs;
                if norm(fc, Inf) < obj.tol
                    break;
                end
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
                %i = i+1;
            end
            %disp(i)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function c = newton(obj, data, ei, mask, rhs, a, b)
            %[f,df] = eval_int_lag_local_newton(obj, data, ei, mask, rhs, cold);
            %c = cold - f./df;
            %i = 0;
            c = 0.5*(a+b);
            rf_flag = true;
            for iter = 1:10
                cold = c;
                [f,df] = eval_int_lag_local_newton(obj, data, ei, mask, rhs, cold);
                step = f./df;
                step(isnan(step)) = 0;
                c  = cold - step;
                I1 = c<a;
                I2 = c>b;
                I3 = ~I1 & ~I2;
                c  = a.*I1 + b.*I2 + c.*I3;
                if ( norm(f, Inf) < obj.tol ) || ( norm(step, Inf) < obj.tol )
                    rf_flag = false;
                    break;
                end
                %i = i+2;
            end
            %disp(i)
            %norm(f, inf)
            if rf_flag
                disp('newton does not converge')
                fc = eval_int_lag_local(obj, data, ei, mask, c) - rhs;
                I1 = (fc < 0);
                I2 = (fc > 0);
                I3 = ~I1 & ~I2;
                a  = I1.*c + I2.*a + I3.*a;
                b  = I1.*b + I2.*c + I3.*b;
                c = regula_falsi(obj, data, ei, mask, rhs, a, b);
            end
        end
        
        %{
        function c = newton(obj, data, ei, mask, rhs, a, b)    
            cold = 0.5*(a+b);
            [f,df] = eval_int_lag_local_newton(obj, data, ei, mask, rhs, cold);
            c = cold - f./df;
            %i = 2;
            while ( norm(c-cold, Inf) >= obj.tol )
                %if ( norm(c-cold, Inf) < obj.tol )
                %    break
                %end
                cold = c;
                [f,df] = eval_int_lag_local_newton(obj, data, ei, mask, rhs, cold);
                c = cold - f./df;
                if norm(f, Inf) < obj.tol
                    break;
                end
                %i = i+2;
            end
            %disp(i)
            %norm(f, inf)
        end
        %}
    end
end