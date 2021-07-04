function def = lag2cheby(type, domain, n_nodes)
%
%Define data structure that maps lagrange polynomials to Chebyshev polynomials
%with preserved boundary values. This function is mainly useful for the squared
%mode--in which the order of the polynomial is doubled. For the non-squared 
%verion, this can be further simplified by directly building the Vandermonde 
%matrix on the original Lagrange nodes. 
%
%%%%%
%Inputs:
%
%type:
%  'Chebyshev1st' or 'Chebyshev2nd', default is 'Chebyshev2nd'
%
%domain:
%  domain of the input polynomial
%
%n_nodes: 
%  number of nodes used in the transformation
%
%
%%%%%
%Output:
%
%def:
%  A data structure contains:
% 
%  nodes:       nodes used for the transformation from Lagrange to Chebyshev
%               we need to evaluate the pdf function on these nodes
%  num_nodes:   number of nodes used in the transformation
%  domain:      domain of the transformation
%  basis2node:  inverse of the vandermonde matrix, transform nodal values
%               to coefficient of 
%  node2basis:  Vandermonde matrix
%
%Tiangang Cui, August, 2019

order = n_nodes - 1;

if n_nodes < 3
    error('Must use more than three nodes')
end

switch type
    case{'Chebyshev1st'}
        tmp     = setup_orthogonal('Chebyshev1st', n_nodes-3);
        int_pts = reference2domain(tmp.ref_nodes, domain);
        def.nodes       = [domain(1); int_pts(:); domain(2)];
        def.basis2node  = eval_orthogonal_basis('Chebyshev1st', domain, order, def.nodes);
        
    case{'Chebyshev2nd'}
        tmp     = setup_orthogonal('Chebyshev2nd', n_nodes-3);
        int_pts = reference2domain(tmp.ref_nodes, domain);
        def.nodes       = [domain(1); int_pts(:); domain(2)];
        def.basis2node  = eval_orthogonal_basis('Chebyshev2nd', domain, order, def.nodes);
end

[L,U] = lu(vpa(def.basis2node));
def.node2basis  = double(U\(L\eye(n_nodes)));
def.num_nodes   = n_nodes;
def.domain      = domain;
def.type        = type;

end