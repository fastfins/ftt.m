function J = eval_rt_jac(firt, r, z)
%Using IRT to sample from the target pdf approximated by firt
%
%%%%
%Inputs:
%
%firt:
%  A data structure built by the function build_irt for evalating the IRT
%
%z:
%  uniform random variables, d x n
%
%%%%%
%Outputs:
%
%r:
%  random variable drawn form the pdf defined by firt
%
%f:
%  pdf function values of those random variables
%
%%%%
%Example:
%
%d = 5;
%a = 0.8;
%
%data = setup_ou_process(d, a);
%func  = @(x) eval_ou_process(data, x); % the pdf
%poly  = setup_oned(30, 'domain', [-5,5]);
%
%opt   = ftt_options('ng_flag', true, 'oned_ref', poly);
%ftt   = build_ftt(func, d, [], opt);
%firt  = build_irt(ftt);
%[r,f] = eval_irt(firt, rand(d, 1E4));
%
%
%Tiangang Cui, August, 2019


d = length(firt.cores);
n = size(r,2);
if size(r,1) ~= d
    error('Input variable must be d x n')
end
J = zeros(d,n*d);

if firt.dir > 0
    % function value of each block
    block_ftt = cell(d,1);
    block_mar = cell(d,1);
    T = cell(d,1);
    block_ftt_d = cell(d,1);
    for k = 1:d
        block_ftt{k} = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(k,:));
        block_mar{k} = eval_oned_core_213(firt.oneds{k}, firt.ys{k}, r(k,:));
        %
        rkm  = size(firt.cores{k}, 1);
        T{k} = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
        %
        block_ftt_d{k} = eval_oned_core_213_deri(firt.oneds{k}, firt.cores{k}, r(k,:));
    end
    %
    F = cell(d,1); % accumulated ftt
    G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
    F{1} = block_ftt{1};
    G{1} = block_mar{1};
    for k = 2:d
        rkm = size(firt.cores{k}, 1);
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
            J(1,ind) = Fm{1}/firt.z;
        else
            J(j,ind) = Fm{j}./Fm{j-1};
        end
        %
        % derivative of the FTT
        drl = block_ftt_d{j};
        % derivative of the FTT, 2nd term, for the d(j+1)/dj term
        mrl = eval_oned_core_213_deri(firt.oneds{j}, firt.ys{j}, r(j,:));
        if j > 1
            rjm = size(firt.cores{j}, 1);
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
            rkm = size(firt.cores{k}, 1);
            nc  = firt.oned_cdfs{k}.num_nodes;
            pk  = reshape(sum( reshape((F{k-1}*T{k}).*(drl*T{k}), n*nc, []), 2), n, nc)';
            % the first term
            J(k,ind) = J(k,ind) + reshape(eval_oned_cdf_deri(firt.oned_cdfs{k}, pk, r(k,:)), 1, []);
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
    
else
    % function value of each block
    block_ftt = cell(d,1);
    block_mar = cell(d,1);
    T = cell(d,1);
    block_ftt_d = cell(d,1);
    for k = 1:d
        block_ftt{k} = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
        block_mar{k} = eval_oned_core_231(firt.oneds{k}, firt.ys{k}, r(k,:));
        rk   = size(firt.cores{k}, 3);
        T{k} = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
        block_ftt_d{k} = eval_oned_core_231_deri(firt.oneds{k}, firt.cores{k}, r(k,:));
    end
    %
    F = cell(d,1); % accumulated ftt
    G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
    F{d} = block_ftt{d}';
    G{d} = block_mar{d}';
    for k = d-1:-1:1
        rk  = size(firt.cores{k}, 3);
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
            J(d,ind) = Fm{d}/firt.z;
        else
            J(j,ind) = Fm{j}./Fm{j+1};
        end
        %
        % derivative of the FTT
        drg = block_ftt_d{j}';
        % derivative of the FTT, 2nd term, for the d(j+1)/dj term
        mrg = eval_oned_core_231_deri(firt.oneds{j}, firt.ys{j}, r(j,:))';
        if j < d
            rj  = size(firt.cores{j}, 3);
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
            rk  = size(firt.cores{k}, 3);
            nc  = firt.oned_cdfs{k}.num_nodes;
            pk  = reshape(sum(reshape((T{k}*F{k+1}).*(T{k}*drg), [], nc*n), 1), nc, n);
            % the first term
            J(k,ind) = J(k,ind) + reshape(eval_oned_cdf_deri(firt.oned_cdfs{k}, pk, r(k,:)), 1, []);
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
