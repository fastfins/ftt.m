function def = setup_lagrange_cdf(oned, local)
%
%Define cdf basis for lagrange polynomials. The key is the global2local matrix
%which maps cdf at global nodes to local nodes. Not proving further details
%here. 
%
%Tiangang Cui, August, 2019

def.num_elems   = oned.num_elems;
def.elem_left   = oned.elem_left;
def.elem_right  = oned.elem_right;
def.elem_sizes  = oned.elem_sizes;
def.gs          = oned.gs;
def.bc          = oned.bc;
def.jac         = oned.jac;
def.domain      = oned.domain;
def.grid_pts    = [def.elem_left def.elem_right(end)];

% setup global nodes
def.num_nodes   = def.num_elems*(local.num_nodes-1)+1;
def.nodes       = zeros(1, def.num_nodes);
for i = 1:def.num_elems
    ind = ( 1:local.num_nodes ) + (local.num_nodes-1)*(i-1);
    def.nodes(ind) = local.nodes*def.elem_sizes(i) + def.elem_left(i);
end

% map the function value y to each local element
if local.num_nodes > 2
    j = local.num_nodes:(local.num_nodes-1):def.num_nodes;
    def.global2local = reshape([reshape(1:(def.num_nodes-1), local.num_nodes-1, def.num_elems); j], [], 1);
else
    def.global2local = reshape([1:(def.num_nodes-1); 2:def.num_nodes], [], 1);
end

def.nodes = def.nodes(:);

end