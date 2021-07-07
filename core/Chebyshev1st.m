classdef Chebyshev1st < spectral
    
    properties
        n(1,:) % only for Chebyshev
    end
    
    methods
        function obj = Chebyshev1st(order, varargin)
            [obj.order,obj.domain] = spectral.process_input(order,varargin{:});
            %
            n = obj.order + 1;
            obj.ref_nodes = reshape( sort( double( cos( vpa(pi)*(2*(1:n)-1)/(2*n) ) ), 'ascend'), [], 1);
            obj.weights = double( ones(size(obj.ref_nodes))*vpa(pi)/n );
            %
            obj.n = reshape(0:obj.order,1,[]);
            obj.normalising = reshape([1,sqrt(2)*ones(1,obj.order)]/sqrt(pi), 1, []);
            %
            obj = post_construction(obj);
        end
        
        
        function [f,w] = eval_ref_basis(obj, x)
            %
            % Evaluate Chebyshev polynomials of the first kind,
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
            % T_n(x) = cos( n * theta )
            % w(x) = 1 / sqrt(1 - x^2)
            %
            theta = real(acos(x(:)));
            f = cos(theta.*obj.n) .* obj.normalising;
            
            if nargout > 1
                w = 1 ./ sqrt( 1 - x(:).^2 );
            end
        end
        
        
        function [f,w] = eval_ref_basis_deri(obj, x)
            %
            % Evaluate derivative of Chebyshev polynomials of the first kind,
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
            % T_n(x) = cos( n * theta )
            % w(x) = 1 / sqrt(1 - x^2)
            % d T_n(x) / dx = n * sin( n * theta) / sin(theta)
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
            f = (sin(theta*obj.n).*obj.n)./sin(theta);
            f = f.*obj.normalising;
            
            if nargout > 1
                w = 1 ./ sqrt( 1 - x(:).^2 );
            end
        end
    end
end