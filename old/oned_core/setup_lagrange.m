function lag = setup_lagrange(grid_pts, local, bc, ghost_size)
%
%Define h-p finite element Lagrange basis, here the order of the basis 
%functions are fixed -- explicitly given by the data structure 'local'. 
%It currently supports no flux and zero essential boundary conditions 
%through the use of ghost elements. 
% 
%In the future, should also consider more complicated boundary type.
% 
%%%%%
%Inputs: 
%
%grid_pts: 
%  A list of points defining the elements
%
%local:    
%  The local Lagrange basis
%
%bc:       
%  'Neumann' or 'Dirichlet'
%
%ghost_size: 
%  Size of the ghost element for handling boundary condition
%
%%%%%
%Output:   
%  lag, a data structure contains:
%
%  nodes:       all the global interpolation points, nx1
%  num_nodes:   number of global nodes
%  num_elems:   number of elements
%  elem_left:   left boudnary of each element, excluding the ghost element,
%               1xm
%  elem_right:  right boudnary of each element, excluding the ghost element
%               1xm
%  elem_sizes:  size of each element, 1xm
%  bc:          type of boundary condition
%  gs:          size of the ghost element
%  mass:        global mass matrix, nxn, sparse
%  mass_L:      Cholesky of the mass matrix, nxn, sparse
%  weights:     weighting factors (integration of each global basis 
%               function), nx1
%  global2local: 
%               this matrix, local.num_nodes x num_elems, contains indices
%               that maps functions evaluated at global nodes to the nodes
%               of local element, a key component used for evaluating
%               Lagrage interpolation
%
%Tiangang Cui, August, 2019


lag.num_elems   = length(grid_pts)-1;
lag.elem_left   = grid_pts(1:end-1);
lag.elem_right  = grid_pts(2:end);
lag.elem_sizes  = (lag.elem_right - lag.elem_left);

% setup global nodes
lag.num_nodes   = lag.num_elems*(local.num_nodes-1)+1;
lag.nodes       = zeros(1, lag.num_nodes);
for i = 1:lag.num_elems
    ind = ( 1:local.num_nodes ) + (local.num_nodes-1)*(i-1);
    lag.nodes(ind) = local.nodes*lag.elem_sizes(i) + lag.elem_left(i);
end

% map the function value y to each local element
if local.num_nodes > 2
    j = local.num_nodes:(local.num_nodes-1):lag.num_nodes;
    lag.global2local = [reshape(1:(lag.num_nodes-1), local.num_nodes-1, lag.num_elems); j]';
else
    lag.global2local = [1:(lag.num_nodes-1); 2:lag.num_nodes]';
end

% setup the weights, mass matrix and its inverse
lag.nodes   = lag.nodes(:);
lag.mass    = zeros(lag.num_nodes);
lag.weights = zeros(lag.num_nodes,1);
lag.jac     = zeros(lag.num_elems,1);
for i = 1:lag.num_elems
    ind = ( 1:local.num_nodes ) + (local.num_nodes-1)*(i-1);
    lag.jac(i) = lag.elem_sizes(i)/(local.domain(2) - local.domain(1));
    lag.mass(ind,ind) = lag.mass(ind,ind) + local.mass*lag.jac(i);
    lag.weights(ind)  = lag.weights(ind) + local.weights(:)*lag.jac(i);
end

% modify the boundary layer
switch bc
    case{'Neumann'}
        tmp_m = ghost_size;
        tmp_w = ghost_size;
    case{'Dirichlet'}
        tmp_m = ghost_size/3;
        tmp_w = ghost_size/2;
    otherwise
        tmp_m = 0;
        tmp_w = 0;
end
lag.mass(1,1)       = lag.mass(1,1) + tmp_m;
lag.mass(end,end)   = lag.mass(end,end) + tmp_m;
lag.weights(1)      = lag.weights(1) + tmp_w;
lag.weights(end)    = lag.weights(end) + tmp_w;

lag.mass    = sparse(0.5*(lag.mass+lag.mass'));
lag.mass_L  = chol(lag.mass)';
%lag.inv_mass = inv(lag.mass);
%lag.inv_mass = 0.5*(lag.inv_mass+lag.inv_mass');

% setup the intergration of each basis

end

