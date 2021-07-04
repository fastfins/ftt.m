function [r,f] = eval_cond_irt(obj, x, z)
% Using SIRT to draw conditional samples from the target pdf approximated
% by FTT
%
%   z   - m x n uniform random variables for generating conditional samples
%
%   x   - given random variable: (d-m) x n
%
%   r   - random variable drawn form the pdf p(.|x) defined by firt
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


d  = length(obj.cores);
nr = size(z,2);
dr = size(z,1);
nx = size(x,2);
dx = size(x,1);

if dx == 0 || dr == 0
    error('dimension of x or the dimension of z should be nonzero')
end

if dr + dx ~= d
    error('dimension of x and the dimension of z mismatch the dimension of the joint')
end

if nx ~= nr
    if nx ~= 1
        error('number of x and the number of z mismatch')
    end
end

r = zeros(size(z));

% first compute the marginal
if obj.dir > 0
    %from the last dimension, marginalise to the first
    %order of the samples (x, r)
    %first evaluate the marginal for x
    frl = ones(nr,1);
    for k = 1:(dx-1)
        rkm = size(obj.cores{k}, 1);
        % evaluate the updated basis function
        T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, x(k,:));
        jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
    k   = dx;
    rkm = size(obj.cores{k}, 1);
    % evaluate the updated basis function
    % for the marginal
    T2  = eval_oned_core_213(obj.oneds{k}, obj.ys{k}, x(k,:));
    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
    ii  = repmat((1:nx)', 1, rkm);
    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
    frl_m   = B*T2;
    fm  = sum(frl_m.^2, 2)';
    
    %for the frl, evaluate the updated basis function
    T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, x(k,:));
    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
    ii  = repmat((1:nx)', 1, rkm);
    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
    frl = B*T2;
    %then generate the conditional samples
    for j = 1:dr
        k = dx+j;
        rkm = size(obj.cores{k}, 1);
        nc  = obj.oned_cdfs{k}.num_nodes;
        
        %bk  = compute_basis(obj.oneds{k}, obj.oned_cdfs{k});
        T1  = reshape(eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), rkm, []);
        pk  = reshape(sum(reshape(frl*T1, nr*nc, []).^2, 2), nr, nc)';
        %
        r(j,:)  = sample_oned_cdf(obj.oned_cdfs{k}, pk, z(j,:));
        %
        % evaluate the updated basis function
        T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, r(j,:));
        jj  = reshape(reshape(1:rkm*nr, rkm, nr)', [], 1);
        ii  = repmat((1:nr)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), nr, rkm*nr);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
    f   = frl'.^2;
else
    %from the first dimension, marginalise to the last
    %order of the samples (r, x)
    %first evaluate the marginal for x
    frg = ones(1,nr);
    for j = dx:-1:2
        k   = dr + j;
        rk  = size(obj.cores{k}, 3);
        % evaluate the updated basis function
        T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, x(j,:));
        ii  = reshape(1:rk*nx, [], 1);
        jj  = reshape(repmat(1:nx, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*nx, nx);
        frg = T2'*B;
    end
    k   = dr+1;
    rk  = size(obj.cores{k}, 3);
    % evaluate the updated basis function
    
    % for marginal
    T2  = eval_oned_core_231(obj.oneds{k}, obj.ys{k}, x(1,:));
    ii  = reshape(1:rk*nx, [], 1);
    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
    B   = sparse(ii, jj, frg(:), rk*nx, nx);
    frg_m   = T2'*B;
    fm  = sum(frg_m.^2, 1);
    
    % for frg
    T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, x(1,:));
    ii  = reshape(1:rk*nx, [], 1);
    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
    B   = sparse(ii, jj, frg(:), rk*nx, nx);
    frg = T2'*B;
    %then generate the conditional samples
    for k = dr:-1:1
        rk  = size(obj.cores{k}, 3);
        nc  = obj.oned_cdfs{k}.num_nodes;
        T1  = reshape(eval_oned_core_213(obj.oneds{k}, obj.ys{k}, obj.oned_cdfs{k}.nodes(:)), [], rk);
        pk  = reshape(sum(reshape(T1*frg, [], nc*nr).^2, 1), nc, nr);
        %
        r(k,:)  = sample_oned_cdf(obj.oned_cdfs{k}, pk, z(k,:));
        %
        % evaluate the updated basis function
        T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, r(k,:));
        ii  = reshape(1:rk*nr, [], 1);
        jj  = reshape(repmat(1:nr, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*nr, nr);
        frg = T2'*B;
    end
    f   = frg.^2;
end

f = f./fm;

end
