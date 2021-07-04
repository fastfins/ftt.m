function z = eval_rt(firt, r)
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
z = zeros(d,n);

if firt.dir > 0
    frl = ones(n,1);
    for k = 1:d
        rkm = size(firt.cores{k}, 1);
        nc  = firt.oned_cdfs{k}.num_nodes;
        
        if firt.ng_flag
            T1  = reshape( eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
            pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)';
        else
            pk  = firt.ys{k}'*frl';
        end
        
        z(k,:)  = eval_oned_cdf(firt.oned_cdfs{k}, pk, r(k,:));
        
        % evaluate the updated basis function
        T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(k,:));
        jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
        ii  = repmat((1:n)', 1, rkm);
        B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
        frl = B*T2;
        %frl(isnan(frl)) = 0;
    end
        
else
    frg = ones(1,n);
    for k = d:-1:1
        rk  = size(firt.cores{k}, 3);
        nc  = firt.oned_cdfs{k}.num_nodes;
        
        if firt.ng_flag
            T1  = reshape( eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
            pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
        else
            pk  = firt.ys{k}*frg;
        end
        
        z(k,:)  = eval_oned_cdf(firt.oned_cdfs{k}, pk, r(k,:));
        
        % evaluate the updated basis function
        T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
        ii  = reshape(1:rk*n, [], 1);
        jj  = reshape(repmat(1:n, rk, 1), [], 1);
        B   = sparse(ii, jj, frg(:), rk*n, n);
        frg = T2'*B;
    end
    
end

end
