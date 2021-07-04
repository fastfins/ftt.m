function T = eval_oned_core_213_deri(oned, core, x)

rkm = size(core, 1);
nn  = oned.num_nodes;
nx  = length(x);

% evaluate the updated basis function
tmp = eval_oned_nodes2x_deri(oned, x(:), reshape(permute(core, [2,1,3]), nn, []));
T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        
end
        
        
        