function J = eval_rt_jac(obj, r, z)
% Evaluate the jacobian of the squared Rosenblatt transport Z = R(X), where 
% Z is the uniform random variable and X is the target random variable. 
%   J = EVAL_RT_JAC(irt, X, Z)
%
%   X - random variable drawn form the pdf defined by SIRT
%   Z - uniform random variables, d x n
%   J - Jacobianx, (d x d) x n, each d x d block is the Jabocian for X(:,j)

d = length(obj.cores);
n = size(r,2);
if size(r,1) ~= d
    error('Input variable must be d x n')
end
J = zeros(d,n*d);

if obj.int_dir > 0 % from left to right
    % function value of each block
    block_ftt = cell(d,1);
    block_mar = cell(d,1);
    T = cell(d,1);
    block_ftt_d = cell(d,1);
    for k = 1:d
        block_ftt{k} = SIRT.eval_oned_core_213(obj.oneds{k}, obj.cores{k}, r(k,:));
        block_mar{k} = SIRT.eval_oned_core_213(obj.oneds{k}, obj.ys{k}, r(k,:));
        %
        rkm  = size(obj.cores{k}, 1);
        T{k} = reshape(SIRT.eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        %
        block_ftt_d{k} = SIRT.eval_oned_core_213_deri(obj.oneds{k}, obj.cores{k}, r(k,:));
    end
    %
    F = cell(d,1); % accumulated ftt
    G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
    F{1} = block_ftt{1};
    G{1} = block_mar{1};
    for k = 2:d
        rkm = size(obj.cores{k}, 1);
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), F{k-1}(:), n, rkm*n);
        F{k} = B*block_ftt{k};
        G{k} = B*block_mar{k};
    end
    %
    Fm = cell(d,1); % accumulated ftt
    for k = 1:d
        Fm{k} = sum(G{k}.^2, 2)';
    end
    %
    % j is the index of the coordinate of differentiation
    for j = 1:d
        ind = (1:d:n*d) + (j-1);
        %
        % the diagonal
        if j == 1
            J(1,ind) = Fm{1}/obj.z;
        else
            J(j,ind) = Fm{j}./Fm{j-1};
        end
        %
        % derivative of the FTT
        drl = block_ftt_d{j};
        % derivative of the FTT, 2nd term, for the d(j+1)/dj term
        mrl = SIRT.eval_oned_core_213_deri(obj.oneds{j}, obj.ys{j}, r(j,:));
        if j > 1
            rjm = size(obj.cores{j}, 1);
            jj  = reshape(reshape(1:rjm*n, rjm, n)', [], 1);
            ii  = repmat((1:n)', 1, rjm);
            B   = sparse(ii(:), jj(:), F{j-1}(:), n, rjm*n);
            drl = B*drl;
            mrl = B*mrl;
        end
        if j < d
            J(j+1,ind) = J(j+1,ind) - sum(G{j}.*mrl,2)'.*z(j+1,:);
        end
        %
        for k = (j+1):d
            % accumulate the j-th block and ealuate the integral
            rkm = size(obj.cores{k}, 1);
            nc  = obj.oned_cdfs{k}.num_nodes;
            pk  = reshape(sum( reshape((F{k-1}*T{k}).*(drl*T{k}), n*nc, []), 2), n, nc)';
            % the first term
            J(k,ind) = J(k,ind) + reshape(eval_cdf_deri(obj.oned_cdfs{k}, pk, r(k,:)), 1, []);
            %
            if k < d
                jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                ii  = repmat((1:n)', 1, rkm);
                B   = sparse(ii(:), jj(:), drl(:), n, rkm*n);
                % the second term, for the d(k+1)/dj term
                mrl = B*block_mar{k};
                J(k+1,ind) = J(k+1,ind) - sum(G{k}.*mrl,2)'.*z(k+1,:);
                % acumulate
                drl = B*block_ftt{k};
            end
            J(k,ind) = 2*J(k,ind)./Fm{k-1};
        end
        %
    end
    
else % from right to left
    % function value of each block
    block_ftt = cell(d,1);
    block_mar = cell(d,1);
    T = cell(d,1);
    block_ftt_d = cell(d,1);
    for k = 1:d
        block_ftt{k} = SIRT.eval_oned_core_231(obj.oneds{k}, obj.cores{k}, r(k,:));
        block_mar{k} = SIRT.eval_oned_core_231(obj.oneds{k}, obj.ys{k}, r(k,:));
        rk   = size(obj.cores{k}, 3);
        T{k} = reshape(SIRT.eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        block_ftt_d{k} = SIRT.eval_oned_core_231_deri(obj.oneds{k}, obj.cores{k}, r(k,:));
    end
    %
    F = cell(d,1); % accumulated ftt
    G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
    F{d} = block_ftt{d}';
    G{d} = block_mar{d}';
    for k = d-1:-1:1
        rk  = size(obj.cores{k}, 3);
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, F{k+1}(:), rk*n, n);
        F{k} = block_ftt{k}'*B;
        G{k} = block_mar{k}'*B;
    end
    %
    Fm = cell(d,1); % accumulated ftt
    for k = 1:d
        Fm{k} = sum(G{k}.^2, 1);
    end
    %
    % j is the index of the coordinate of differentiation
    for j = d:-1:1
        ind = (1:d:n*d) + (j-1);
        %
        % the diagonal
        if j == d
            J(d,ind) = Fm{d}/obj.z;
        else
            J(j,ind) = Fm{j}./Fm{j+1};
        end
        %
        % derivative of the FTT
        drg = block_ftt_d{j}';
        % derivative of the FTT, 2nd term, for the d(j+1)/dj term
        mrg = SIRT.eval_oned_core_231_deri(obj.oneds{j}, obj.ys{j}, r(j,:))';
        if j < d
            rj  = size(obj.cores{j}, 3);
            ii  = reshape(1:rj*n, [], 1);
            jj  = reshape(repmat(1:n, rj, 1), [], 1);
            B   = sparse(ii, jj, F{j+1}(:), rj*n, n);
            drg = drg*B;
            mrg = mrg*B;
        end
        if j > 1
            J(j-1,ind) = J(j-1,ind) - sum(G{j}.*mrg,1).*z(j-1,:);
        end
        %
        for k = (j-1):-1:1
            % accumulate the j-th block and ealuate the integral
            rk  = size(obj.cores{k}, 3);
            nc  = obj.oned_cdfs{k}.num_nodes;
            pk  = reshape(sum(reshape((T{k}*F{k+1}).*(T{k}*drg), [], nc*n), 1), nc, n);
            % the first term
            J(k,ind) = J(k,ind) + reshape(eval_cdf_deri(obj.oned_cdfs{k}, pk, r(k,:)), 1, []);
            %
            if k > 1
                ii  = reshape(1:rk*n, [], 1);
                jj  = reshape(repmat(1:n, rk, 1), [], 1);
                B   = sparse(ii, jj, drg(:), rk*n, n);
                % the second term, for the d(k+1)/dj term
                mrg = block_mar{k}'*B;
                J(k-1,ind) = J(k-1,ind) - sum(G{k}.*mrg,1).*z(k-1,:);
                % acumulate
                drg = block_ftt{k}'*B;
            end
            J(k,ind) = 2*J(k,ind)./Fm{k+1};
        end
        %
    end
end
end