function fx = eval_ftt_block(ftt, x, dir)
%Evaluate the ftt function for either the first or last k variables. 
%The output is horizontally aligned.
%
%Inputs: 
%
%ftt:
%  A given function tensor train
%
%x: 
%  input variables, d x n. 
%
%dir:
%  direction of the evaluation      
%
%Outputs:
%
%fx:
%  function values at x, m x n
%
%Tiangang Cui, August, 2019

d   = length(ftt.cores);
k   = size(x, 1);
nx  = size(x, 2);

if dir > 0
    % start from the 1st dimension
    fx  = ones(nx,1);
    
    % all the intermediate dimensions, except the last dimension
    for j = 1:min(k,d)
        nj  = ftt.oneds{j}.num_nodes;
        rjm = size(ftt.cores{j}, 1);
        %
        if j < d || (size(ftt.cores{j}, 3) > 1 && size(ftt.cores{j}, 4) == 1)
            tmp = reshape(permute(ftt.cores{j}, [2,1,3]), nj, []);
        else
            % collapse the third dimension = 1
            tmp = reshape(permute(reshape(ftt.cores{j}, rjm, nj, []), [2,1,3]), nj, []);
        end

        T   = reshape(permute(reshape( eval_oned_nodes2x(ftt.oneds{j},x(j,:),tmp), nx, rjm, []), [2,1,3]), rjm*nx, []);
        
        % how to speed up this part?
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        B   = sparse(ii(:), jj(:), fx(:), nx, rjm*nx);
        
        fx  = B*T;
    end
else
    % start from the last dimension
    xind = k:-1:1;
    tind = d:-1:1;
    fx   = ones(1,nx);
    % all the intermediate dimensions, except the first dimension
    % need to walk through d-k+1 dimensions
    for i = 1:min(k,d)
        j = tind(i);
        
        nj  = ftt.oneds{j}.num_nodes;
        rj  = size(ftt.cores{j}, 3);
        %
        if j > 1 || (size(ftt.cores{j}, 1) > 1 && size(ftt.cores{j}, 4) == 1)
            tmp = reshape(permute(ftt.cores{j}, [2,3,1]), nj, []);
        else
            % collapse the first dimension = 1
            tmp = reshape(ftt.cores{j}, nj, []);
        end
        
        T   = reshape(permute(reshape( eval_oned_nodes2x(ftt.oneds{j},x(xind(i),:),tmp), nx, rj, []), [2,1,3]), rj*nx, []);
        
        % how to speed up this part?
        ii  = reshape(1:rj*nx, [], 1);
        jj  = reshape(repmat(1:nx, rj, 1), [], 1);
        B   = sparse(ii, jj, fx(:), rj*nx, nx);
        
        fx  = T'*B;
    end
end


end