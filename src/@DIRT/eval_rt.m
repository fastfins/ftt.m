function [z, logf] = eval_rt(obj, x, k)
% Evaluate the deep Rosenblatt transport Z = T(X), where Z is the reference
% random variable and X is the target random variable.
%   [z, logf] = EVAL_RT(dirt, x, k)
%
%   x    - random variable drawn form the pdf defined by DIRT
%   k    - number of layers in the evaluations
%   z    - reference random variables, d x n
%   logf - log of the DIRT density

if nargin <= 2
    k = obj.n_layers;
else
    k = min(k, obj.n_layers);
end
z = x;
logf = zeros(1, size(x,2));
for l = 1:k
    % ref density, skip the first one
    if l > 1
        logfz = log_pdf(obj.diag, z);
    else
        logfz = 0;
    end
    % evaluate rt
    u = eval_rt(obj.irts{l}, z);
    % eval pdf (det Jacobian)
    f = eval_pdf(obj.irts{l}, z);
    % evaluate the diagonal transform
    z = eval_icdf(obj.diag, u);
    % update density
    logf = logf + log(f) - logfz;
end

end
