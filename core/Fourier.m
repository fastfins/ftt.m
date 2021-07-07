classdef Fourier < spectral
    
    properties 
        m
        n
        is
        ic
        c
    end
    
    methods
        function obj = Fourier(order, varargin)
            [obj.order,obj.domain] = spectral.process_input(order,varargin{:});
            %
            n = obj.order*2 + 2;
            obj.ref_nodes = reshape( sort( (2/n)*(1:n) - 1, 'ascend'), [], 1);
            obj.weights = ones(size(obj.ref_nodes))*(2/n);
            %
            obj.m   = obj.order+1;
            obj.n   = obj.m*2;
            obj.is  = 2:obj.m;
            obj.ic  = (obj.m+1):(obj.n-1);
            obj.c   = reshape((1:obj.order)*pi, 1, []);
            %
            obj = post_construction(obj);
        end
        
        
        function [f,w] = eval_ref_basis(obj, x)
            %
            tmp = x(:).*obj.c;     
            f = [ones(size(x(:)))*sqrt(0.5), sin(tmp), cos(tmp), ...
                cos( x(:)*(obj.m*pi) )*sqrt(0.5)];
            %{
            f = zeros(length(x), obj.n);
            f(:,1)  = ones(size(x(:)))*sqrt(0.5);
            f(:,obj.is) = sin( tmp );
            f(:,obj.ic) = cos( tmp );
            f(:,obj.n)  = cos( x(:)*(obj.m*pi) )*sqrt(0.5);
            %}
            w = ones(size(x));
        end
        
        function [f, w] = eval_ref_basis_deri(obj, x)
            %
            tmp = x(:).*obj.c;
            f = [zeros(length(x),1), cos(tmp).*obj.c, -sin(tmp).*obj.c, ...
                -sin(x(:)*(obj.m*pi))*sqrt(0.5)*(obj.m*pi)];
            %{
            f = zeros(length(x), obj.n);
            f(:,obj.is) =  cos( tmp ).*obj.c;
            f(:,obj.ic) = -sin( tmp ).*obj.c;
            f(:,obj.n)  = -sin( x(:)*(obj.m*pi) )*sqrt(0.5)*(obj.m*pi);
            %}
            w = ones(size(x));
        end
        
        function b = eval_ref_int_basis(obj, x)
            %
            tmp = x(:).*obj.c;
            b = [x(:)*sqrt(0.5), -cos(tmp)./obj.c, sin(tmp)./obj.c, ...
                sin(x(:)*(obj.m*pi))*sqrt(0.5)/(pi*obj.m)];
            %{
            b = zeros(length(x), obj.n);
            b(:,1)  = x(:)*sqrt(0.5);
            b(:,obj.is) = - cos( tmp )./ obj.c;
            b(:,obj.ic) = sin( tmp )  ./ obj.c;
            b(:,obj.n)  = sin( x(:)*(obj.m*pi) )*sqrt(0.5)/(pi*obj.m);
            %}
        end
        
        function [b,db] = eval_ref_int_basis_newton(obj, x)
            %
            tmp = x(:).*obj.c;
            ct  = cos(tmp);
            st  = sin(tmp);
            %
            b = [x(:)*sqrt(0.5), -ct./obj.c, st./obj.c, ...
                sin(x(:)*(obj.m*pi))*sqrt(0.5)/(pi*obj.m)];
            %
            db = [ones(size(x(:)))*sqrt(0.5), st, ct, ...
                cos(x(:)*(obj.m*pi))*sqrt(0.5)];
        end
    end
end