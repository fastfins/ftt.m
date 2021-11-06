classdef Lagrange1 < piecewise
    
    properties
        local_mass(2,2) = [2, 1; 1, 2]/6
        local_weights(1,2) = [1, 1]/2
        local_domain(1,2) = [0, 1]
    end
    
    methods
        function obj = Lagrange1(num_elems, varargin)
            obj@piecewise(1, num_elems, varargin{:});
            %
            obj.nodes = obj.grid;
            obj.num_nodes = length(obj.grid);
            %
            obj.jac     = obj.elem_size;
            obj.mass    = zeros(obj.num_nodes);
            obj.weights = zeros(obj.num_nodes,1);
            for i = 1:obj.num_elems
                ind = [i, i+1];
                obj.mass(ind,ind) = obj.mass(ind,ind) + obj.local_mass*obj.jac;
                obj.weights(ind)  = obj.weights(ind)  + obj.local_weights(:)*obj.jac;
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
            %
            obj.mass    = sparse(0.5*(obj.mass+obj.mass'));
            obj.mass_R  = chol(obj.mass);
            obj.int_W   = reshape(obj.weights, 1, []);
            %
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %{
        function f = eval(obj, f_at_nodes, x)
            tau = eps; % safe guard thershold
            n   = length(x);
            m   = size(f_at_nodes,2);
            f   = zeros(n,m);
            
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
            if sum(mask_left) > 0
                switch obj.bc
                    case{'Dirichlet'}
                        txl = max(obj.domain(1), x(mask_left));
                        d = left_node-obj.domain(1);
                        f(mask_left,:)  = ((txl-obj.domain(1))/d).*f_at_nodes(1,:);
                    otherwise
                        f(mask_left,:)  = repmat(f_at_nodes(1,:), sum(mask_left), 1);
                end
            end
            if sum(mask_right) > 0
                switch obj.bc
                    case{'Dirichlet'}
                        txr = min(obj.domain(2), x(mask_right));
                        d = right_node-obj.domain(2);
                        f(mask_right,:) = ((txr-obj.domain(2))/d).*f_at_nodes(obj.num_nodes,:);
                    otherwise
                        f(mask_right,:) = repmat(f_at_nodes(obj.num_nodes,:), sum(mask_right), 1);
                end
            end
            
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                ind = ceil((tmp_x-left_node)./obj.elem_size);
                ind(ind==0) = 1;
                %
                % f_left  = f_at_nodes(ind);
                % f_right = f_at_nodes(ind+1);
                local_x = (reshape(tmp_x, [], 1) - reshape(obj.grid(ind), [], 1))./obj.elem_size;
                %
                f(mask_inside,:) = f_at_nodes(ind,:).*(1-local_x) + f_at_nodes(ind+1,:).*local_x;
            end
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function bas = eval_basis(obj, x)
            tau = eps; % safe guard thershold
            n   = length(x);
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
            roi = [find(mask_left); find(mask_right)];
            coi = [ones(sum(mask_left),1); ones(sum(mask_right),1)*obj.num_nodes];
            switch obj.bc
                case{'Dirichlet'}
                    txl = max(obj.domain(1), x(mask_left));
                    txr = min(obj.domain(2), x(mask_right));
                    val = [reshape( ((txl-obj.domain(1))/(left_node-obj.domain(1))), [], 1);...
                        reshape(((txr-obj.domain(2))/(right_node-obj.domain(2))), [], 1)];
                otherwise
                    val = [ones(sum(mask_left),1); ones(sum(mask_right),1)];
            end
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                ind = ceil((tmp_x-left_node)./obj.elem_size);
                ind(ind==0) = 1;
                %
                % map each x into local coordinate
                local_x = (reshape(tmp_x, 1, []) - reshape(obj.grid(ind), 1, []))./obj.elem_size;
                %left:  ind: 1-local_x
                %right: ind+1: local_x
                %
                coi = [coi; ind(:); ind(:)+1];
                roi = [roi; repmat(reshape(find(mask_inside),[],1), 2, 1)];
                val = [val; 1-local_x(:); local_x(:)];
            end
            bas = sparse(roi, coi, val, n, obj.num_nodes);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function bas = eval_basis_deri(obj, x)
            tau = eps; % safe guard thershold
            n   = length(x);
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
                    roi = [find(mask_left); find(mask_right)];
                    coi = [ones(sum(mask_left),1); ones(sum(mask_right),1)*obj.num_nodes];
                    val = [ones(sum(mask_left),1)./(left_node-obj.domain(1));...
                        ones(sum(mask_right),1)./(right_node-obj.domain(2))];
                otherwise
                    roi = [];
                    coi = [];
                    val = [];
            end
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                ind = ceil((tmp_x-left_node)./obj.elem_size);
                ind(ind==0) = 1;
                %
                coi = [coi; ind(:); ind(:)+1];
                roi = [roi; repmat(reshape(find(mask_inside),[],1), 2, 1)];
                val = [val; -ones(length(ind), 1)./obj.elem_size; ones(length(ind),1)./obj.elem_size];
            end
            bas = sparse(roi, coi, val, n, obj.num_nodes);
        end
    end
end