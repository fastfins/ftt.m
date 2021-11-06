function [f, g] = pullback(obj, func, z)
% Evaluate the pullback function $T^\sharp\pi$ for a given deep Rosenblatt
% transport X=T(Z), where Z is the reference random variable and X is the
% target random variable.
%   [f,g] = PULLBACK(dirt, func, z)
%
%   func - function handle of the target function, returns minus log
%          likelihood and minus log prior, and their gradients
%   z    - reference random variables, d x n
%   f    - minus log target function value
%   g    - gradient of minus log target function value

if nargout == 1
    [x,logft] = eval_irt(obj, z);
    % compute the minus log likelihood and minus log prior
    [mllkds,mlps] = func(x);
    % compute the reference density at z
    logfz = log_pdf(obj.diag, z);
    %f = -(mllkds+mlps) - logft + logfz;
    f = (mllkds+mlps) + logft - logfz;
else
    [x,logft,gz,Juz,Jux] = eval_irt(obj, z);
    % compute the minus log likelihood and minus log prior
    [mllkds,mlps,gmllkds,gmlps] = func(x);
    % compute the reference density at z
    [logfz,glfz] = log_pdf(obj.diag, z);
    f = (mllkds+mlps) + logft - logfz;
    g = DIRT.backtrack(Juz, Jux, (gmllkds+gmlps)) + gz - glfz;
end

end