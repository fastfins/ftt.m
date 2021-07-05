function [r,f] = eval_irt(obj, z)
% Evaluate squared IRT r = T(z), where z is uniform
%
%   z   - uniform random variables, d x n
%
%   r   - random variable drawn form the pdf defined by SIRT
%
%   f   - pdf function values of those random variables
%
% Example:
%
% % setup the model problem
% d = 5;
% a = 0.8;
% data = setup_ou_process(d, a);
% func  = @(x) eval_ou_process(data, x); % the pdf
% % setup polynomial
% poly  = Legendre(30, [-5,5]);
% % options for building FTT
% opt   = FTToption('sqrt_flag', true);
% % build FTT
% ftt   = FTT(func, d, poly, opt);
% build IRT
% irt  = SIRT(ftt);
% ealuate
% [r,f] = eval_irt(irt, rand(d, 1E4));
%


d = length(obj.cores);
[dz,n] = size(z);
r = zeros(dz,n);

if obj.marginal_direction > 0 % from left to right
    frl = ones(n,1);
    for k = 1:dz
        rkm = size(obj.cores{k}, 1);
        nc  = obj.oned_cdfs{k}.num_nodes;
        T1  = reshape( SIRT.eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)';
        %
        r(k,:)  = invert_cdf(obj.oned_cdfs{k}, pk, z(k,:));
        %
        % evaluate the updated basis function
        T2  = SIRT.eval_oned_core_213(obj.oneds{k}, obj.cores{k}, r(k,:));
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
    if dz < d
        f   = sum((frl*obj.ms{dz}).^2, 2)'/obj.z;
    else
        f   = frl'.^2/obj.z;
    end
else % from right to left
    frg = ones(1,n);
    ie  = (d-dz)+1;
    for k = d:-1:ie
        ck  = k-ie+1;
        rk  = size(obj.cores{k}, 3);
        nc  = obj.oned_cdfs{k}.num_nodes;
        T1  = reshape( SIRT.eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
        %
        r(ck,:)  = invert_cdf(obj.oned_cdfs{k}, pk, z(ck,:));
        %
        % evaluate the updated basis function
        T2  = SIRT.eval_oned_core_231(obj.oneds{k}, obj.cores{k}, r(ck,:));
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*n, n);
        frg = T2'*B;
    end
    if dz < d
        f  = sum((obj.ms{ie}*frg).^2, 1)/obj.z;
    else
        f   = frg.^2/obj.z;
    end
    
end

end