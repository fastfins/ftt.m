function bas = eval_lagrange_basis_deri(def, x)
%Evaluate basis of piecewise lagrange polynomials at given x coordinates
%except two ghost cells on the left and right, all the other elements have
%same order and same size
%
%Tiangang Cui, August, 2019

tau = eps; % safe guard thershold
n   = length(x);
bas = zeros(n,def.num_nodes);

if sum(def.domain(1)-tau > x | def.domain(2)+tau < x)
    disp('warning: points outside of the domain')
end

% build mask to remove points that are too close to boudnary
left_node   = def.elem_left(1);
right_node  = def.elem_right(def.num_elems);
mask_left   = left_node-eps >= x(:);
mask_right  = right_node <= x(:);
mask_inside = ~(mask_left | mask_right);

% boundary
%{
switch def.bc
    case{'Dirichlet'}
        bas(mask_left,1) = ((x(mask_left)-def.domain(1))/(left_node-def.domain(1)))*1;
        bas(mask_right,def.num_nodes) = ((x(mask_right)-def.domain(2))/(right_node-def.domain(2)))*1;
    otherwise
        bas(mask_left,1) = 1;
        bas(mask_right,def.num_nodes) = 1;
end
%}
switch def.bc
    case{'Dirichlet'}
        bas(mask_left,1) = 1/(left_node-def.domain(1));
        bas(mask_right,def.num_nodes) = 1/(right_node-def.domain(2));
    otherwise
        %bas(mask_left,1) = 0;
        %bas(mask_right,def.num_nodes) = 0;
end

if sum(mask_inside) > 0
    
    tmp_x   = x(mask_inside);
    % find the element indices for each x
    
    ind     = sum(def.elem_left(:)'-eps < tmp_x(:), 2);
    
    % map each x into local coordinate
    local_x = (reshape(tmp_x, [], 1) - reshape(def.elem_left(ind), [], 1))./reshape(def.elem_sizes(ind), [], 1);
    
    % evaluate the barycentric formula
    diff    = repmat(local_x(:), 1, def.local.num_nodes) - repmat(def.local.nodes(:)', length(tmp_x), 1);
    % stablise
    diff(abs(diff)<tau) = tau;
    tmp_m1 = repmat(def.local.omega, length(tmp_x), 1) ./ diff;
    tmp_m2 = repmat(def.local.omega, length(tmp_x), 1) ./ (diff.^2);
    %original function
    %lbs     = tmp_m./sum(tmp_m, 2);
    %
    %
    a = 1./sum(tmp_m1, 2);
    b = sum(tmp_m2, 2).*(a.^2);
    lbs = (tmp_m1.*b - tmp_m2.*a)./reshape(def.jac(ind), [], 1);
    
    % embed lbs into the global grid
    % def.global2local(ind,:) are the col indices
    % repmat(find(mask_inside), 1, def.local.num_nodes)  are the row indices
    %
    
    coi = def.global2local(ind,:);
    roi = repmat(find(mask_inside), 1, def.local.num_nodes);
    ii  = (coi-1)*n + roi;
    
    % evaluation of the internal interpolation
    
    bas(ii(:))  = lbs(:);
end

end



