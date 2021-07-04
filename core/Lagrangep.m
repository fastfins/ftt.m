classdef Lagrangep < piecewise
    
    properties
        local LagrangeRef
        global2local(:,:) 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function obj = Lagrangep(order, num_elems, varargin)
            obj@piecewise(order, num_elems, varargin{:});
            %
            if order == 1
                disp('should use Lagrange1')
            end
            %
            obj.local = LagrangeRef(obj.order+1);
            
            % setup global nodes
            obj.num_nodes   = obj.num_elems*(obj.local.num_nodes-1)+1;
            obj.nodes       = zeros(1, obj.num_nodes);
            for i = 1:obj.num_elems
                ind = ( 1:obj.local.num_nodes ) + (obj.local.num_nodes-1)*(i-1);
                obj.nodes(ind) = obj.local.nodes*obj.elem_size + obj.grid(i);
            end
            
            % map the function value y to each local element
            if obj.local.num_nodes > 2
                j = obj.local.num_nodes:(obj.local.num_nodes-1):obj.num_nodes;
                obj.global2local = [reshape(1:(obj.num_nodes-1), obj.local.num_nodes-1, obj.num_elems); j]';
            else
                obj.global2local = [1:(obj.num_nodes-1); 2:obj.num_nodes]';
            end
            
            % setup the weights, mass matrix and its inverse
            obj.jac     = obj.elem_size/(obj.local.domain(2) - obj.local.domain(1));
            obj.nodes   = obj.nodes(:);
            obj.mass    = zeros(obj.num_nodes);
            obj.weights = zeros(obj.num_nodes,1);
            for i = 1:obj.num_elems
                ind = ( 1:obj.local.num_nodes ) + (obj.local.num_nodes-1)*(i-1);
                obj.mass(ind,ind) = obj.mass(ind,ind) + obj.local.mass*obj.jac;
                obj.weights(ind)  = obj.weights(ind) + obj.local.weights(:)*obj.jac;
            end
            
            % modify the boundary layer
            switch obj.bc
                case{'Neumann'}
                    tmp_m = obj.gs;
                    tmp_w = obj.gs;
                case{'Dirichlet'}
                    tmp_m = obj.gs/3;
                    tmp_w = obj.gs/2;
                otherwise
                    tmp_m = 0;
                    tmp_w = 0;
            end
            obj.mass(1,1)       = obj.mass(1,1) + tmp_m;
            obj.mass(end,end)   = obj.mass(end,end) + tmp_m;
            obj.weights(1)      = obj.weights(1) + tmp_w;
            obj.weights(end)    = obj.weights(end) + tmp_w;
            
            obj.mass    = sparse(0.5*(obj.mass+obj.mass'));
            obj.mass_L  = chol(obj.mass)';
            %
            %
            obj.name = ['Lagrangep'];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval(obj, f_at_nodes, x)
            bas = eval_basis(obj, x(:));
            f = bas*f_at_nodes;
            
            %{
            tau = eps; % safe guard thershold
            n   = length(x);
            f   = zeros(n,1);
            
            if sum(obj.domain(1)-tau > x | obj.domain(2)+tau < x)
                disp('warning: points outside of the domain')
            end
            
            % build mask to remove points that are too close to boudnary
            left_node   = obj.grid(1);
            right_node  = obj.grid(end);
            mask_left   = left_node-eps > x(:);
            mask_right  = right_node+eps < x(:);
            mask_inside = ~(mask_left | mask_right);
            
            % boundary
            switch obj.bc
                case{'Dirichlet'}
                    d = left_node-obj.domain(1);
                    f(mask_left)    = ((x(mask_left)-obj.domain(1))/d)*f_at_nodes(1);
                    d = right_node-obj.domain(2);
                    f(mask_right)   = ((x(mask_right)-obj.domain(2))/d)*f_at_nodes(obj.num_nodes);
                otherwise
                    f(mask_left)    = f_at_nodes(1);
                    f(mask_right)   = f_at_nodes(obj.num_nodes);
            end
            
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                % ind = sum(obj.elem_left(:)'-eps <= tmp_x(:), 2);
                %
                ind = ceil((tmp_x(:)-left_node)./obj.elem_size);
                ind(ind==0) = 1;
                %
                % map the function value y to each local element
                % num_elem x num_local_nodes
                local_f = reshape(f_at_nodes(obj.global2local), size(obj.global2local));
                % local_f(ind,:) gives the local function value for each local element
                
                % map each x into local coordinate
                local_x = (reshape(tmp_x, 1, []) - reshape(obj.grid(ind), 1, []))./obj.elem_size;
                %end
                
                % evaluate the barycentric formula
                diff    = repmat(local_x(:), 1, obj.local.num_nodes) - repmat(obj.local.nodes(:)', length(tmp_x), 1);
                % stablise
                diff(abs(diff)<tau) = tau;
                tmp_m   = repmat(obj.local.omega, length(tmp_x), 1) ./ diff;
                
                % evaluation of the internal interpolation
                f(mask_inside)  = sum(local_f(ind,:).*tmp_m, 2)./sum(tmp_m, 2);
            end
            %}
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [bas, w] = eval_basis(obj, x)
            
            tau = eps; % safe guard thershold
            n   = length(x);
            bas = zeros(n,obj.num_nodes);
            
            if sum(obj.domain(1)-tau > x | obj.domain(2)+tau < x)
                disp('warning: points outside of the domain')
            end
            
            % build mask to remove points that are too close to boudnary
            left_node   = obj.grid(1);
            right_node  = obj.grid(end);
            mask_left   = left_node-eps >= x(:);
            mask_right  = right_node <= x(:);
            mask_inside = ~(mask_left | mask_right);
            
            % boundary
            switch obj.bc
                case{'Dirichlet'}
                    txl = max(obj.domain(1), x(mask_left));
                    txr = min(obj.domain(2), x(mask_right));
                    bas(mask_left,1) = ((txl-obj.domain(1))/(left_node-obj.domain(1)))*1;
                    bas(mask_right,obj.num_nodes) = ((txr-obj.domain(2))/(right_node-obj.domain(2)))*1;
                otherwise
                    bas(mask_left,1) = 1;
                    bas(mask_right,obj.num_nodes) = 1;
            end
            
            
            if sum(mask_inside) > 0
                
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                
                ind = ceil((tmp_x(:)-left_node)./obj.elem_size);
                ind(ind==0) = 1;
                
                % map each x into local coordinate
                local_x = (reshape(tmp_x, 1, []) - reshape(obj.grid(ind), 1, []))./obj.elem_size;
                
                % evaluate the barycentric formula
                diff    = repmat(local_x(:), 1, obj.local.num_nodes) - repmat(obj.local.nodes(:)', length(tmp_x), 1);
                % stablise
                diff(abs(diff)<tau) = tau;
                tmp_m   = repmat(obj.local.omega, length(tmp_x), 1) ./ diff;
                lbs     = tmp_m./sum(tmp_m, 2);
                
                % embed lbs into the global grid
                % obj.global2local(ind,:) are the col indices
                % repmat(find(mask_inside), 1, obj.local.num_nodes)  are the row indices
                %
                
                coi = obj.global2local(ind,:);
                roi = repmat(find(mask_inside), 1, obj.local.num_nodes);
                ii  = (coi-1)*n + roi;
                
                % evaluation of the internal interpolation
                
                bas(ii(:))  = lbs(:);
            end
            
            w = ones(size(x(:)));
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [bas, w] = eval_basis_deri(obj, x)
            
            tau = eps; % safe guard thershold
            n   = length(x);
            bas = zeros(n,obj.num_nodes);
            
            if sum(obj.domain(1)-tau > x | obj.domain(2)+tau < x)
                disp('warning: points outside of the domain')
            end
            
            % build mask to remove points that are too close to boudnary
            left_node   = obj.grid(1);
            right_node  = obj.grid(end);
            mask_left   = left_node-eps >= x(:);
            mask_right  = right_node <= x(:);
            mask_inside = ~(mask_left | mask_right);
            
            switch obj.bc
                case{'Dirichlet'}
                    bas(mask_left,1) = 1/(left_node-obj.domain(1));
                    bas(mask_right,obj.num_nodes) = 1/(right_node-obj.domain(2));
                otherwise
                    %bas(mask_left,1) = 0;
                    %bas(mask_right,obj.num_nodes) = 0;
            end
            
            if sum(mask_inside) > 0
                
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                
                ind = ceil((tmp_x(:)-left_node)./obj.elem_size);
                ind(ind==0) = 1;
                
                % map each x into local coordinate
                local_x = (reshape(tmp_x, 1, []) - reshape(obj.grid(ind), 1, []))./obj.elem_size;
                
                % evaluate the barycentric formula
                diff    = repmat(local_x(:), 1, obj.local.num_nodes) - repmat(obj.local.nodes(:)', length(tmp_x), 1);
                % stablise
                diff(abs(diff)<tau) = tau;
                tmp_m1 = repmat(obj.local.omega, length(tmp_x), 1) ./ diff;
                tmp_m2 = repmat(obj.local.omega, length(tmp_x), 1) ./ (diff.^2);
                %original function
                %lbs     = tmp_m./sum(tmp_m, 2);
                %
                %
                a = 1./sum(tmp_m1, 2);
                b = sum(tmp_m2, 2).*(a.^2);
                lbs = (tmp_m1.*b - tmp_m2.*a)./reshape(obj.jac(ind), [], 1);
                
                % embed lbs into the global grid
                % obj.global2local(ind,:) are the col indices
                % repmat(find(mask_inside), 1, obj.local.num_nodes)  are the row indices
                %
                
                coi = obj.global2local(ind,:);
                roi = repmat(find(mask_inside), 1, obj.local.num_nodes);
                ii  = (coi-1)*n + roi;
                
                % evaluation of the internal interpolation
                
                bas(ii(:))  = lbs(:);
            end
            
            w = ones(size(x(:)));
            
        end
        
    end
end