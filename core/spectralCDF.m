classdef spectralCDF < onedCDF
    % spectralCDF class 
    %
    % For Fourier basis, FourierCDF is used. For other spectral polynomials
    % in bounded domains, we first transform the polynomial to the 2nd 
    % Chebyshev basis, and then apply the inversion. See Chebyshev2ndCDF.
    %
    % Before applying root findings, a grid search based on sampling_nodes
    % is applied to locate the left and right boundary of root finding. 
    %
    % See also ChebyshevCDF and FourierCDF.
    
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
        
        function [f,df] = eval_int_newton(obj, coef, c, x)
            [x,J] = domain2reference(obj,x(:));
            [b,db] = eval_ref_int_basis_newton(obj, x(:));
            b = b.*J;
            %
            if size(coef,2) > 1
                if size(coef,2) == length(x)
                    f  = reshape(sum(b .*coef',2) - c, [], 1);
                    df = reshape(sum(db.*coef',2), [], 1);
                else
                    error('Error: dimenion mismatch')
                end
            else
                f  = reshape((b *coef), [], 1) - c;
                df = reshape((db*coef), [], 1);
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
                %r(mask2) = regula_falsi(obj, coef, cdf_base+rhs(mask2), a(:), b(:));
                r(mask2) = newton(obj, coef, cdf_base+rhs(mask2), a(:), b(:));
                if sum( r(mask2)>b(:) | r(mask2)<a(:) ) ~=0
                    warning('newton failed');
                    r(mask2) = regula_falsi(obj, coef, cdf_base+rhs(mask2), a(:), b(:));
                end
            else
                %r(mask2) = regula_falsi(obj, coef(:,mask2), ...
                %    reshape(cdf_base(mask2),[],1)+rhs(mask2), a(:), b(:));
                r(mask2) = newton(obj, coef(:,mask2), reshape(cdf_base(mask2),[],1)+rhs(mask2), a(:), b(:));
                if sum( r(mask2)>b(:) | r(mask2)<a(:) ) ~=0
                    warning('newton failed');
                    r(mask2) = regula_falsi(obj, coef(:,mask2), ...
                        reshape(cdf_base(mask2),[],1)+rhs(mask2), a(:), b(:));
                end
            end
            %
            r = reshape(r, size(xi));
            r(isnan(r)) = 0.5*(obj.domain(1) + obj.domain(2));
            r(isinf(r)) = 0.5*(obj.domain(1) + obj.domain(2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function c = regula_falsi(obj, coef, rhs, a, b)
            fa = eval_int(obj,coef,a)-rhs;
            fb = eval_int(obj,coef,b)-rhs;
            if sum(sign(fb.*fa) ~= -1)
                disp('Root finding: initial guesses on one side')
            end
            c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
            cold = inf;
            %i = 2;
            while ( norm(c-cold, Inf) > obj.tol )
                cold = c;
                fc  = eval_int(obj,coef,c)-rhs;
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
        
        function c = newton(obj, coef, rhs, a, b)
            %i = 0;
            c = 0.5*(a+b);
            rf_flag = true;
            for iter = 1:10
                cold = c;
                [f,df] = eval_int_newton(obj, coef, rhs, cold);
                step = f./df;
                step(isnan(step)) = 0;
                c = cold - step;
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
                fc = eval_int(obj,coef,c)-rhs;
                I1 = (fc < 0);
                I2 = (fc > 0);
                I3 = ~I1 & ~I2;
                a  = I1.*c + I2.*a + I3.*a;
                b  = I1.*b + I2.*c + I3.*b;
                c = regula_falsi(obj, data, ei, mask, rhs, a, b);
            end
        end
        
        %{
        function c = newton(obj, coef, rhs, a, b)
            cold = 0.5*(a+b);
            [f,df] = eval_int_newton(obj, coef, rhs, cold);
            c = cold - f./df;
            %i = 2;
            while ( norm(c-cold, Inf) >= obj.tol )
                cold = c;
                [f,df] = eval_int_newton(obj, coef, rhs, cold);
                if norm(f, Inf) < obj.tol
                    break;
                end
                c = cold - f./df;
                %norm(f, inf)
                %i = i+2;
            end
            %disp(i)
        end
        %}
    end
    
end