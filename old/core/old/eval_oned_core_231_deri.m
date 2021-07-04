function T = eval_oned_core_231_deri(oned, core, x)

rk  = size(core, 3);
nn  = oned.num_nodes;
nx  = length(x);

% evaluate the updated basis function
tmp = eval_oned_nodes2x_deri(oned, x(:), reshape(permute(core, [2,3,1]), nn, []));
T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        
end
      
