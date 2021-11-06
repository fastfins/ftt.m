classdef Chebyshev2nd < spectral
    
    properties
        n(1,:) % only for Chebyshev
    end
    
    methods
        function obj = Chebyshev2nd(order, varargin)
            [obj.order,obj.domain] = spectral.process_input(order,varargin{:});
            %
            n = obj.order + 1;
            obj.ref_nodes = reshape( sort( double( cos( vpa(pi)*(1:n)/(n+1) ) ), 'ascend'), [], 1);
            obj.weights = double( sin( vpa(pi)*(1:n)/(n+1)).^2*vpa(pi)/(n+1) );
            %
            obj.n = reshape(0:obj.order,1,[]);
            obj.normalising = sqrt(2/pi);  
            %
            obj = post_construction(obj);
        end
        
        function [f,w] = eval_ref_basis(obj, x)
            %
            % Evaluate Chebyshev polynomials of the 2nd kind,
            % for all input x, up to order n
            %
            % Inputs:
            % x:    n_pts
            %
            % Output:
            % f:    function outputs for each order at each x, n_pts x (n+1)
            % w:    weight function at each x, n_pts x 1
            %
            % x = cos(theta), or theta = acos(x), x in [-1, 1]
            %
            % U_n(theta) = cos( (n+1) * theta )  / cos(theta)
            % w(x) = sqrt(1 - x^2)
            %
            theta = real(acos(x(:)));
            % deal with end points
            f = sin(theta.*(obj.n+1)) ./ (sin(theta)./obj.normalising);
            
            mask = abs(x+1) < eps;
            if sum(mask) > 0
                f(mask,:) = repmat(((obj.n+1).*(-1).^obj.n).*obj.normalising, sum(mask), 1);
            end
            
            mask = abs(x-1) < eps;
            if sum(mask) > 0
                f(mask,:) = repmat((obj.n+1).*obj.normalising, sum(mask), 1);
            end
            
            if nargout > 1
                w = sqrt( 1 - x(:).^2 );
            end
        end
        
        function [f,w] = eval_ref_basis_deri(obj, x)
            %
            % Evaluate derivative of Chebyshev polynomials of the second kind,
            % for all input x, up to order n
            %
            % Inputs:
            % x:    n_pts
            %
            % Output:
            % f:    function outputs for each order at each x, n_pts x (n+1)
            % w:    weight function at each x, n_pts x 1
            %
            % x = cos(theta), or theta = acos(x), x in [-1, 1]
            % ( d theta / dx ) = inv(dx/dtheta) = inv( - sin(theta))
            %
            % sin(theta)^2 = (1-x^2)
            %
            % U_n(theta) = sin( (n+1) * theta )  / sin(theta)
            % w(x) = sqrt(1 - x^2)
            % d U_n(x) / dx = - (n+1) * cos((n+1)*theta) / sin(theta)^2 + sin((n+1)*theta) cos(theta) / sin(theta)^3
            %               = ( - (n+1) * cos((n+1)*theta) + x sin((n+1)*theta) / sin(theta) ) / (1-x^2)
            %               = ( (n+1) * cos((n+1)*theta) - x sin((n+1)*theta) / sin(theta) ) / (x^2-1)
            %
            theta = real(acos(x(:)));
            
            % deal with end points
            mask = abs(x+1) < eps;
            if sum(mask) > 0
                theta(mask) = pi;
            end
            mask = abs(x-1) < eps;
            if sum(mask) > 0
                theta(mask) = 0;
            end
            
            f = ( cos(theta*(obj.n+1)).*(obj.n+1) - sin(theta*(obj.n+1)).*(x(:)./sin(theta)) )./(x(:).^2-1);
            f = f.*obj.normalising;
            
            if nargout > 1
                w = sqrt( 1 - x(:).^2 );
            end
        end
        
        function b = eval_ref_int_basis(obj, x)
            %
            theta = real(acos(x));
            % the normalising constant used
            % int(U_n) = T_(n+1) / (n+1), n = 0, ..., order
            b = cos( theta.*(obj.n+1) ) .* (obj.normalising./(obj.n+1));
        end
        
        function [b,db] = eval_ref_int_basis_newton(obj, x)
            %
            % Evaluate Chebyshev polynomials of the 2nd kind,
            % for all input x, up to order n
            %
            % Inputs:
            % x:    n_pts
            %
            % Output:
            % f:    function outputs for each order at each x, n_pts x (n+1)
            % w:    weight function at each x, n_pts x 1
            %
            % x = cos(theta), or theta = acos(x), x in [-1, 1]
            %
            % U_n(theta) = cos( (n+1) * theta )  / cos(theta)
            % w(x) = sqrt(1 - x^2)
            %
            theta = real(acos(x(:)));
            b  = cos(theta.*(obj.n+1)) .* (obj.normalising./(obj.n+1));
            db = sin(theta.*(obj.n+1)) ./ (sin(theta)./obj.normalising);
            
            % deal with end points
            
            mask = abs(x+1) < eps;
            if sum(mask) > 0
                db(mask,:) = repmat(((obj.n+1).*(-1).^obj.n).*obj.normalising, sum(mask), 1);
            end
            
            mask = abs(x-1) < eps;
            if sum(mask) > 0
                db(mask,:) = repmat((obj.n+1).*obj.normalising, sum(mask), 1);
            end
        end
    end
    
end