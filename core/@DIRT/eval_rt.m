function [z, logf] = eval_rt(obj, x)
% Evaluate the deep Rosenblatt transport Z = T(X), where Z is the reference
% random variable and X is the target random variable.
%   [Z, logf] = EVAL_RT(dirt, X)
%
%   X    - random variable drawn form the pdf defined by DIRT
%   Z    - reference random variables, d x n
%   logf - log of the DIRT density

z = x;
logf = zeros(1, size(x,2));
for l = 1:obj.n_layers
    % ref density, skip the first one
    if l > 1
        logfz = eval_log_pdf(obj.diag, z);
    else
        logfz = 0;
    end
    % evaluate rt
    u = eval_rt(obj.irts{l}, z);
    % eval pdf (det Jacobian)
    f = eval_pdf(obj, z);
    % evaluate the diagonal transform
    z = eval_icdf(obj.diag, u);
    % update density
    logf = logf + log(f) - logfz;
end

end
