function f = eval_lagrange(def, f_at_nodes, x)
%Evaluate piecewise lagrange polynomials at given x coordinates
%except two ghost cells on the left and right, all the other elements have
%same order and same size.
%
%Tiangang Cui, August, 2019

tau = eps; % safe guard thershold
n   = length(x);
f   = zeros(n,1);

if sum(def.domain(1)-tau > x | def.domain(2)+tau < x)
    disp('warning: points outside of the domain')
end

% build mask to remove points that are too close to boudnary
left_node   = def.elem_left(1);
right_node  = def.elem_right(def.num_elems);
mask_left   = left_node-eps > x(:);
mask_right  = right_node+eps < x(:);
mask_inside = ~(mask_left | mask_right);

% boundary
switch def.bc
    case{'Dirichlet'}
        d = left_node-def.domain(1);
        f(mask_left)    = ((x(mask_left)-def.domain(1))/d)*f_at_nodes(1);
        d = right_node-def.domain(2);
        f(mask_right)   = ((x(mask_right)-def.domain(2))/d)*f_at_nodes(def.num_nodes);
    otherwise
        f(mask_left)    = f_at_nodes(1);
        f(mask_right)   = f_at_nodes(def.num_nodes);
end

if sum(mask_inside) > 0
    tmp_x   = x(mask_inside);
    % find the element indices for each x
    ind     = sum(def.elem_left(:)'-eps <= tmp_x(:), 2);
    % map the function value y to each local element
    % num_elem x num_local_nodes
    local_f = reshape(f_at_nodes(def.global2local), size(def.global2local));
    % local_f(ind,:) gives the local function value for each local element
    
    % map each x into local coordinate
    %if def.num_elems == 1
    %    local_x = (reshape(tmp_x, 1, []) - def.elem_left)/def.elem_sizes;
    %else
    local_x = (reshape(tmp_x, 1, []) - reshape(def.elem_left(ind), 1, []))./reshape(def.elem_sizes(ind), 1,[]);
    %end
    
    % evaluate the barycentric formula
    diff    = repmat(local_x(:), 1, def.local.num_nodes) - repmat(def.local.nodes(:)', length(tmp_x), 1);
    % stablise
    diff(abs(diff)<tau) = tau;
    tmp_m   = repmat(def.local.omega, length(tmp_x), 1) ./ diff;
    
    % evaluation of the internal interpolation
    f(mask_inside)  = sum(local_f(ind,:).*tmp_m, 2)./sum(tmp_m, 2);
end

end



