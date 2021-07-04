classdef LagrangeRef
    %
    %Define the reference Lagrange basis, in the reference domain [0,1]
    %This function should not be explicitly used
    %
    %%%%%
    %Input:
    %
    %type:
    %  Type of the interpolation points, options are 'Jacobi' and 'Chebyshev' (2nd)
    %
    %n:
    %  Number of interpolation points, should be greater than or equal to 2
    %
    %%%%%
    %Output:
    %  lag, a data structure contains:
    %
    %  domain:     the reference domain
    %  nodes:      all the interpolation points, 1xn
    %  num_nodes:  number of nodes
    %  omega:      barycentric weights of the Lagrange polynomial
    %  mass:       reference mass matrix
    %  weights:    reference weighting factors (integration of each basis
    %              function), 1xn
    %
    %Tiangang Cui, August, 2019
    
    properties
        domain(1,2) = [0, 1]
        order 
        nodes(1,:) 
        omega(1,:) 
        num_nodes 
        mass(:,:) 
        weights(1,:) 
    end
    
    methods
        function obj = LagrangeRef(n)
            if n < 2
                error('We need more than two points to define Lagrange interpolation')
            end
            
            obj.nodes       = zeros(1,n);
            obj.nodes(1)    = 0;
            obj.nodes(end)  = 1;
            obj.num_nodes   = n;
            
            if n > 2
                order = n-3;
                Jacob = Jacobi11(order, [-1,1]);
                % to the interval [0, 1]
                obj.nodes(2:n-1) = 0.5*(Jacob.ref_nodes+1);
            end
            
            % compute the local omega coefficients
            obj.omega = zeros(1,n);
            for j = 1:n
                ind     = true(n,1);
                ind(j)  = false;
                obj.omega(j) = 1./prod( obj.nodes(j) - obj.nodes(ind) );
            end
            
            % define the mass matrix
            I = eye(n);
            obj.mass = zeros(n);
            for i = 1:n
                for j = 1:n
                    fij = @(x) eval(obj, I(:,i), x).*eval(obj, I(:,j), x);
                    obj.mass(i,j) = integral(fij, obj.domain(1), obj.domain(2));
                end
            end
            % setup the intergration of each basis
            obj.weights = zeros(1,n);
            for i = 1:n
                fi = @(x) eval(obj, I(:,i), x);
                obj.weights(i) = integral(fi, obj.domain(1), obj.domain(2));
            end
            
        end
        
        function f = eval(obj, f_at_x, x)
            tau = eps;
            m   = length(x);
            n   = obj.num_nodes;
            f   = zeros(m,1);
            %
            outside = obj.nodes(1)-tau >= x | obj.nodes(n)+tau <= x;
            inside  = ~outside;
            if sum(outside)
                disp('warning: points outside of the domain')
                f(outside) = 0;
            end
            %
            tmp_x   = x(inside);
            diff    = repmat(tmp_x(:), 1, n) - repmat(obj.nodes(:)', length(tmp_x), 1);
            %
            % stablise
            diff(abs(diff)<tau) = tau;
            tmp_m   = repmat(obj.omega(:)', length(tmp_x), 1) ./ diff;
            %
            % evaluation of the internal interpolation
            f(inside)   = sum(repmat(f_at_x(:)',length(tmp_x),1).*tmp_m, 2)./sum(tmp_m, 2);
            %
            f = reshape(f, size(x));
        end
        
    end
end