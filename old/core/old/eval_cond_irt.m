function [r,f] = eval_cond_irt(firt, x, z)
%Using IRT to draw conditional samples from the target pdf approximated by 
%firt
%
%%%%
%Inputs: 
%
%firt:
%  A data structure built by the function build_irt for evalating the IRT
%
%z: 
%  uniform random variables, m x n, used for generating conditional samples
%
%x: 
%  given random variable: (d-m) x n
%
%%%%%
%Outputs:
%
%r: 
%  random variable drawn form the pdf p(.|x) defined by firt
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


d  = length(firt.cores);
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
if firt.dir > 0
    %from the last dimension, marginalise to the first
    %order of the samples (x, r)
    %first evaluate the marginal for x
    
    frl = ones(nr,1);
    for k = 1:(dx-1)
        rkm = size(firt.cores{k}, 1);
        % evaluate the updated basis function
        T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, x(k,:));
        jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
    
    k   = dx;
    rkm = size(firt.cores{k}, 1);
    % evaluate the updated basis function
    % for the marginal
    if firt.ng_flag
        T2  = eval_oned_core_213(firt.oneds{k}, firt.ys{k}, x(k,:));
        jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
        frl_m   = B*T2;
        fm  = sum(frl_m.^2, 2)';
    else
        fxk = eval_oned_nodes2x(firt.oneds{k}, x(k,:), firt.ys{k}');
        % mutiply together
        fm  = sum(frl.*fxk, 2)';
    end
    
    %for the frl, evaluate the updated basis function
    T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, x(k,:));
    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
    ii  = repmat((1:nx)', 1, rkm);
    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
    frl = B*T2;
    
    %then generate the conditional samples
    for j = 1:dr
        k = dx+j;
        rkm = size(firt.cores{k}, 1);
        nc  = firt.oned_cdfs{k}.num_nodes;
        
        %bk  = compute_basis(firt.oneds{k}, firt.oned_cdfs{k});
        if firt.ng_flag
            T1  = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
            pk  = reshape(sum(reshape(frl*T1, nr*nc, []).^2, 2), nr, nc)';
        else
            pk  = firt.ys{k}'*frl';
        end
        
        r(j,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(j,:));
        
        % evaluate the updated basis function
        T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(j,:));
        jj  = reshape(reshape(1:rkm*nr, rkm, nr)', [], 1);
        ii  = repmat((1:nr)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), nr, rkm*nr);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end

    if firt.ng_flag
        f   = frl'.^2;
    else
        f   = frl';
    end
else
    %from the first dimension, marginalise to the last
    %order of the samples (r, x)
    %first evaluate the marginal for x
    
    frg = ones(1,nr);
    for j = dx:-1:2
        k   = dr + j;
        rk  = size(firt.cores{k}, 3);
        % evaluate the updated basis function
        T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, x(j,:));
        ii  = reshape(1:rk*nx, [], 1);
        jj  = reshape(repmat(1:nx, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*nx, nx);
        frg = T2'*B;
    end
    
    k   = dr+1;
    rk  = size(firt.cores{k}, 3);
    % evaluate the updated basis function
    
    % for marginal
    if firt.ng_flag
        T2  = eval_oned_core_231(firt.oneds{k}, firt.ys{k}, x(1,:));
        ii  = reshape(1:rk*nx, [], 1);
        jj  = reshape(repmat(1:nx, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*nx, nx);
        frg_m   = T2'*B;
        fm  = sum(frg_m.^2, 1);
    else
        fxk = eval_oned_nodes2x(firt.oneds{k}, x(1,:), firt.ys{k});
        % mutiply together
        fm  = sum(fxk'.*frg, 1);
    end
    
    % for frg
    T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, x(1,:));
    ii  = reshape(1:rk*nx, [], 1);
    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
    B   = sparse(ii, jj, frg(:), rk*nx, nx);
    frg = T2'*B;
    
    %then generate the conditional samples
    for k = dr:-1:1
        rk  = size(firt.cores{k}, 3);
        nc  = firt.oned_cdfs{k}.num_nodes;
        
        if firt.ng_flag
            T1  = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
            pk  = reshape(sum(reshape(T1*frg, [], nc*nr).^2, 1), nc, nr);
        else
            pk  = firt.ys{k}*frg;
        end
        
        r(k,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(k,:));
        
        % evaluate the updated basis function
        T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
        ii  = reshape(1:rk*nr, [], 1);
        jj  = reshape(repmat(1:nr, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*nr, nr);
        frg = T2'*B;
    end
    
    if firt.ng_flag
        f   = frg.^2;
    else
        f   = frg;
    end
end

f = f./fm;

end
