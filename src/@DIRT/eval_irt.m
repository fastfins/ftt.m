function [x,logf,gz,Juz,Jux] = eval_irt(obj, z, k)
% Evaluate the deep Rosenblatt transport X = T(Z), where Z is the reference
% random variable and X is the target random variable.
%   [x,logf,Juz,Jux] = EVAL_IRT(dirt, z, k)
%
%   z    - reference random variables, d x n
%   k    - number of layers in the evaluations
%   x    - random variable drawn form the pdf defined by DIRT
%   logf - log of the DIRT density
%   gz   - gradient of log of the DIRT density
%   Juz  - cell array of the Jacobian of the diagonal map, U = CDF(Z), at
%          each layer, dimension of each: d x n
%   Jux  - cell array of the Jacobian of the SRT (CDF map), U = R(X), at
%          each layer, dimension of each: d x (dn)

if nargin <= 2
    k = obj.n_layers;
else
    k = min(k, obj.n_layers);
end
x = z;
logf = zeros(1, size(z,2));

if nargout <= 2
    for l = k:-1:1
        % evaluate the diagonal transform
        u = eval_cdf(obj.diag, x);
        % evaluate sirt
        [x, ft] = eval_irt(obj.irts{l}, u);
        % evaluate the denominator of the Jacobian
        if l > 1
            logfx = log_pdf(obj.diag, x);
        else
            logfx = 0;
        end
        % update density
        logf = logf + log(ft) - logfx;
    end
else
    Juz = cell(k,1);
    Jux = cell(k,1);
    gx  = cell(k,1);
    for l = k:-1:1
        % evaluate the diagonal transform
        [u,Juz{l}] = eval_cdf(obj.diag, x);
        % evaluate sirt
        [x,ft,glt] = eval_irt(obj.irts{l}, u);
        % compuate Jacobian
        Jux{l} = eval_rt_jac(obj.irts{l}, x, u);
        % evaluate the denominator of the Jacobian
        if l > 1
            [logfx,glfx] = log_pdf(obj.diag, x);
        else
            logfx = 0;
            glfx = 0;
        end
        gx{l} = glt - glfx;
        % update density
        logf = logf + log(ft) - logfx;
    end
    
    [d,n] = size(z);
    gz = zeros(d,n);
    for l = 1:k
        gz = gz + gx{l};
        for i = 1:n
            ind = (i-1)*d + (1:d);
            gz(:,i) = Juz{l}(:,i) .* ( tril(Jux{l}(:,ind))' \ gz(:,i) );
        end
    end

end

end