function [x,logf,Juz,Jux] = eval_irt(obj, z)
% Evaluate the deep Rosenblatt transport X = T(Z), where Z is the reference 
% random variable and X is the target random variable. 
%   [x,logf,Juz,Jux] = EVAL_IRT(dirt, z)
%
%   z    - reference random variables, d x n
%   x    - random variable drawn form the pdf defined by DIRT
%   logf - log of the DIRT density
%   Juz  - cell array of the Jacobian of the diagonal map, U = CDF(Z), at 
%          each layer, dimension of each: d x n
%   Jux  - cell array of the Jacobian of the SRT (CDF map), U = R(X), at 
%          each layer, dimension of each: d x (dn)

x = z;
logf = zeros(1, size(z,2));
if nargout <= 2
    for l = obj.n_layers:-1:1
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
    Juz = cell(obj.n_layers,1);
    Jux = cell(obj.n_layers,1);
    for l = obj.n_layers:-1:1
        % evaluate the diagonal transform
        [u, Juz{l}] = eval_cdf(obj.diag, x);
        % evaluate sirt
        [x, ft] = eval_irt(obj.irts{l}, u);
        % compuate Jacobian
        Jux{l} = eval_rt_jac(obj.irts{l}, x, u);
        % evaluate the denominator of the Jacobian
        if l > 1
            logfx = log_pdf(obj.diag, x);
        else
            logfx = 0;
        end
        % update density
        logf = logf + log(ft) - logfx;
    end
end

end