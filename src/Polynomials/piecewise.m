classdef piecewise < oned
% piecewise class   - Superclass for Lagrange1 (piecewise linear basis) and 
%                     Larangep (piecewise high order Lagrange)
%
% Constructors:
%   * Lagrange1(n, domain, 'ghost_size', gs, 'bc', type)
%   * Lagrangep(order, n, domain, 'ghost_size', gs, 'bc', type)
%
%   order       - (Required) order of the polynomial 
%   n           - (Required) number of elements
%   domain      - (Optional) default is [0,1]
%   ghost_size  - The parameter 'ghost_size' specify a boundary layer with
%                 size gs
%   bc          - The parameter 'bc' specify the type of the boundary layer
%                 to either 'Dirichlet' or 'Neumann'. With 'Dirichlet', the
%                 boudnary layer linearly decays to zero. With 'Neumann',
%                 the boundary layer takes a constant value.
%
% See also PIECEWISECDF
    
    properties
        gs % ghostsize
        bc % boundary condition
        grid(:,1) 
        num_elems 
        elem_size 
        jac 
        mass(:,:) 
        weights(:,1) 
    end
    
    methods (Static)
        function [order,num_elems,domain,bc,gs] = process_input(order, num_elems, varargin)
            defaultDomain = [0, 1];
            defaultBC   = 'Dirichlet';
            expectedBC  = {'Dirichlet','Neumann','none'};
            defaultGhostSize = 0;
            %
            p = inputParser;
            % valid order
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
            addRequired(p,'order',validScalarPosNum);
            % number of lagrange elements
            addRequired(p,'num_elems',@(x)isnumeric(x) && isscalar(x) && (x >= 1));
            % validate domain
            validDomain = @(x) length(x)==2 && x(2)>x(1);
            addOptional(p,'domain',defaultDomain,validDomain);
            % lagrange boundary condition
            addParameter(p,'bc',defaultBC,@(x)any(validatestring(x,expectedBC)));
            % lagrange boundary ghost cell size
            addParameter(p,'ghost_size',defaultGhostSize,validScalarPosNum);
            %
            p.KeepUnmatched = true;
            parse(p,order,num_elems,varargin{:});
            %
            order = p.Results.order;
            num_elems = p.Results.num_elems;
            domain = p.Results.domain;
            bc = p.Results.bc;
            gs = p.Results.ghost_size;            
            %
            if gs < eps
                bc  = 'none';
            end
            if (domain(1) + gs*2) > domain(2)
                error('ghost cell too large')
            end
        end
    end
    
    methods
        function obj = piecewise(order, num_elems, varargin)
            [obj.order,obj.num_elems,obj.domain,obj.bc,obj.gs] = ...
                piecewise.process_input(order, num_elems, varargin{:});
            obj.grid = linspace(obj.domain(1)+obj.gs, obj.domain(2)-obj.gs, obj.num_elems+1);
            obj.elem_size = (obj.grid(end)-obj.grid(1))/double(obj.num_elems);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function x = sample_domain(obj, n)
            x = rand(1,n)*(obj.grid(end)-obj.grid(1))+obj.grid(1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [ind,B,interp_atx] = point_selection(obj, int_method, H)
            % Build the cross indices
            rold = size(H,1)/obj.num_nodes;
            %
            switch int_method
                case {'QDEIM'}
                    [~,~,e] = qr(H', 'vector');
                    ind = e(1:size(H,2));
                    interp_atx = H(ind,:);
                    B =  H/interp_atx;
                case {'DEIM'}
                    disp('Not implemented')
                case {'MaxVol'}
                    [ind,B] = maxvol(H);
                    interp_atx = H(ind,:);
            end
            if cond(interp_atx) > 1E5
                disp('warning: poor condition number in the interpolation')
            end
        end
        
    end
    
end
