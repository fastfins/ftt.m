function logf = log_pdf(obj, x, k)
% Evaluate the deep Rosenblatt transport Z = T(X), where Z is the reference
% random variable and X is the target random variable.
%   logf = LOG_PDF(dirt, x, k)
%
%   x    - random variable drawn form the pdf defined by DIRT
%   k    - number of layers in the evaluations
%   logf - log of the DIRT density

if nargin <= 2
    k = obj.n_layers;
else
    k = min(k, obj.n_layers);
end
[~,logf] = eval_rt(obj, x, k);

end
