classdef Fourier < spectral
    
    methods
        function obj = Fourier(order, varargin)
            [obj.order,obj.domain] = spectral.process_input(order,varargin{:});
            %
            obj.name = 'Fourier';
            %
            n = obj.order*2 + 2;
            obj.ref_nodes = reshape( sort( (2/n)*(1:n) - 1, 'ascend'), [], 1);
            obj.weights = ones(size(obj.ref_nodes))*(2/n);
            %
            obj.name = 'Fourier';
            %
            obj = post_construction(obj);
        end
        
        
        function [f,w] = eval_ref_basis(obj, x)
            %
            m   = obj.order+1;
            n   = m*2;
            is  = 2:m;
            ic  = (m+1):(n-1);
            
            f = zeros(length(x), n);
            f(:,1)  = ones(size(x(:)))*sqrt(0.5);
            f(:,is) = sin( x(:)*(1:obj.order)*pi );
            f(:,ic) = cos( x(:)*(1:obj.order)*pi );
            f(:,n)  = cos( x(:)*(m*pi) )*sqrt(0.5);
            
            w = ones(size(x));
        end
        
        function [f, w] = eval_ref_basis_deri(obj, x)
            %
            m   = obj.order+1;
            n   = m*2;
            is  = 2:m;
            ic  = (m+1):(n-1);
            
            f = zeros(length(x), n);
            c = reshape((1:obj.order)*pi, 1, []);
            f(:,is) =  cos( x(:)*c ).*c;
            f(:,ic) = -sin( x(:)*c ).*c;
            f(:,n)  = -sin( x(:)*(m*pi) )*sqrt(0.5)*(m*pi);
            
            w = ones(size(x));
        end
        
        function b = eval_ref_int_basis(obj, x)
            %
            m   = obj.order+1;
            n   = m*2;
            is  = 2:m;
            ic  = (m+1):(n-1);
            %
            b = zeros(length(x), n);
            b(:,1)  = x(:)*sqrt(0.5);
            b(:,is) = - cos( x(:)*(1:obj.order)*pi )./ (pi*(1:obj.order));
            b(:,ic) = sin( x(:)*(1:obj.order)*pi )  ./ (pi*(1:obj.order));
            b(:,n)  = sin( x(:)*(m*pi) )*sqrt(0.5)/(pi*m);
            %
        end
    end
end