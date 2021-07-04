function F = eval_lagrange_cdf_local(oned_cdf, data, ei, mask, x_mask)
%%Evaluate the cdf function defined by the cdf data strucuture. This only
%deals with a single element in the h-p FEM discretisation of pdf. 
%
%Inputs:
%
%oned_cdf:
%  reference oned cdf basis
%
%data:
%  data used for holding the transformation from pdf to cdf
%
%ei:
%  indicate which local element the inversion of cdf should take place
%
%mask:
%  a mask matrix filter out all the indices falls outside of the local
%  element
%
%x_mask:
%  masked x values
%
%Outputs:
%
%F:       cdf values for x_mask
%
%Tiangang Cui, August, 2019

domains = [reshape(oned_cdf.elem_left(ei(mask)),[],1), reshape(oned_cdf.elem_right(ei(mask)),[],1)];
if data.size == 1
    b   = eval_oned_int_basis(oned_cdf.cheby, domains, oned_cdf.order, x_mask);
    tmp = reshape(sum(b.*data.coef(:,ei(mask))',2), size(x_mask));

    F   = (tmp - reshape(data.base(ei(mask)), size(tmp)))./data.norm + reshape(data.cdf_grid(ei(mask)), size(tmp));
else
    ii  = reshape(1:data.size, size(ei));
    j1  = ii(mask) + (ei(mask)-1)*data.size;
    coe = reshape(data.coef, oned_cdf.local.num_nodes, []);
    
    b   = eval_oned_int_basis(oned_cdf.cheby, domains, oned_cdf.order, x_mask);
    tmp = reshape(sum(b.*coe(:,j1)',2), size(x_mask));
    
    j2  = (ii(mask)-1)* oned_cdf.num_elems    + ei(mask);
    j3  = (ii(mask)-1)*(oned_cdf.num_elems+1) + ei(mask);
    F   = (tmp-reshape(data.base(j2),size(tmp)))./reshape(data.norm(mask),size(tmp)) + reshape(data.cdf_grid(j3),size(tmp));
end

end
