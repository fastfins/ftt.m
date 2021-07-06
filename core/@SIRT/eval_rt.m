function z = eval_rt(obj, r)
% Evaluate the squared Rosenblatt transport Z = R(X), where Z is the uniform 
% random variable and X is the target random variable. 
%   Z = EVAL_RT(irt, X)
%
%   X - random variable drawn form the pdf defined by SIRT
%   Z - uniform random variables, d x n


d = length(obj.cores);
[dr,n] = size(r);
z = zeros(dr,n);

if obj.int_dir > 0 % from left to right
    frl = ones(n,1);
    for k = 1:dr
        rkm = size(obj.cores{k}, 1);
        nc  = obj.oned_cdfs{k}.num_nodes;
        %
        T1  = reshape( SIRT.eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)';
        %
        z(k,:)  = eval_cdf(obj.oned_cdfs{k}, pk, r(k,:));
        %
        % evaluate the updated basis function
        T2  = SIRT.eval_oned_core_213(obj.oneds{k}, obj.cores{k}, r(k,:));
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
else % from right to left
    frg = ones(1,n);
    ie  = (d-dr)+1;
    for k = d:-1:ie
        ck  = k-ie+1;
        rk  = size(obj.cores{k}, 3);
        nc  = obj.oned_cdfs{k}.num_nodes;
        T1  = reshape( SIRT.eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
        %
        z(ck,:)  = eval_cdf(obj.oned_cdfs{k}, pk, r(ck,:));
        %
        % evaluate the updated basis function
        T2  = SIRT.eval_oned_core_231(obj.oneds{k}, obj.cores{k}, r(ck,:));
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*n, n);
        frg = T2'*B;
    end
end

end