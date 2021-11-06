classdef spectral < oned
% spectral class    - Superclass for Fourier, Chebyshev1st, Chebyshev2nd,
%                     and recurr subclasses. The recurr class implements
%                     Legendre, Jabobi11, Hermite, and Laguerre subclasses
%
% Constructors:
%   * poly(order, domain)
%
%   poly        - Choose one from Fourier, Chebyshev1st, Chebyshev2nd, 
%                 Legendre, Jabobi11, Hermite, and Laguerre. 
%   order       - (Required) order of the polynomial 
%   domain      - (Optional) default is [0,1]
%
% See also SPECTRALCDF
    
    properties
        ref_nodes(:,1) 
        weights(:,1) 
        normalising(1,:) 
        omegas(:,1) 
        basis2node(:,:) 
        node2basis(:,:) 
    end
    
    methods (Abstract)
        eval_ref_basis(obj)
        eval_ref_basis_deri(obj)
    end
    
    methods (Static)
        function [order,domain] = process_input(order, varargin)
            defaultDomain = [0, 1];
            p = inputParser;
            % valid order
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
            addRequired(p,'order',validScalarPosNum);
            % validate domain
            validDomain = @(x) (length(x)==2) && (x(2)>x(1));
            addOptional(p,'domain',defaultDomain,validDomain);
            %
            parse(p,order,varargin{:});
            order = p.Results.order;
            domain = p.Results.domain;
        end
    end
    
    methods        
        function obj = post_construction(obj)
            obj.num_nodes = length(obj.ref_nodes);
            [obj.nodes,J] = reference2domain(obj, obj.ref_nodes);
            [obj.basis2node,obj.omegas] = eval_ref_basis(obj, obj.ref_nodes);
            obj.omegas  = reshape(obj.omegas, size(obj.ref_nodes));
            obj.omegas  = obj.omegas/J;
            obj.weights = reshape(obj.weights, size(obj.ref_nodes));
            obj.node2basis  = obj.basis2node'*diag(obj.weights);
            %
            obj.mass_R  = sqrt(obj.weights./obj.omegas).*obj.basis2node;
            obj.int_W   = reshape(obj.weights(:)./obj.omegas(:), 1, [])*obj.basis2node;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x,x2z] = reference2domain(obj, z)
            x2z = 0.5*(obj.domain(:,2)-obj.domain(:,1));
            mid = 0.5*(obj.domain(:,2)+obj.domain(:,1));
            x   = z(:).*x2z+mid;
            x   = reshape(x,size(z));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [z,x2z] = domain2reference(obj, x)
            x2z = 0.5*(obj.domain(:,2)-obj.domain(:,1));
            mid = 0.5*(obj.domain(:,2)+obj.domain(:,1));
            z   = (x(:)-mid)./x2z;
            z   = reshape(z,size(x));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function x = sample_domain(obj, n)
            x = rand(1,n)*(obj.domain(2)-obj.domain(1))+obj.domain(1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [f,w] = eval_basis(obj, x)
            [x,J] = domain2reference(obj, x(:));
            [f,w] = eval_ref_basis(obj, x(:));
            w = w./J;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [f,w] = eval_basis_deri(obj, x)
            [x,J] = domain2reference(obj, x(:));
            [f,w] = eval_ref_basis_deri(obj, x(:));
            w = w./J;
            f = f./J;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [ind,B,interp_atx] = point_selection(obj, int_method, H)
            % Build the cross indices
            %
            rold  = size(H,1)/obj.num_nodes;
            nodes = reshape(obj.basis2node*reshape(H, obj.num_nodes, []), obj.num_nodes, rold, []);
            nodes = reshape(nodes, rold*obj.num_nodes, []);
            switch int_method
                case {'QDEIM'}
                    [~,~,e] = qr(nodes', 'vector');
                    ind = e(1:size(nodes,2));
                    interp_atx = nodes(ind,:);
                    B =  H/interp_atx;
                case {'MaxVol'}
                    [ind,~] = maxvol(nodes);
                    interp_atx = nodes(ind,:);
                    B =  H/interp_atx;
            end
            if cond(interp_atx) > 1E5
                disp('warning: poor condition number in the interpolation')
            end
        end
    end
    
end
