function fx = eval_ftt(ftt, x, varargin)
%Evaluate the ftt function. The output is horizontally aligned.
%
%Inputs: 
%
%ftt:
%  A given function tensor train
%
%x: 
%  input variables, d x n. 
%
%Optional input:
%
%pt:
%  permuation index, where ftt is built for f(y), y is the permutation of x
%  s.t. x = y(p) and y = x(pt)
%
%  for a given p, pt can be obtained by pt(p) = 1:length(p);
%      
%
%Outputs:
%
%fx:
%  function values at x, m x n
%
%Tiangang Cui, August, 2019

% permute the input
if length(varargin) > 1
    error('too many optional arguments')
elseif length(varargin) == 1
    pt = varargin{1};
    x = x(pt,:);
end

fx  = eval_ftt_block(ftt, x, ftt.direction);

if ftt.direction > 0
    fx  = fx';
end

if ftt.ng_flag
    fx  = fx.^2;
end

end
