function [d,rs,ns] = ftt_size(ftt)
%Find the ranks of the ftt cores and the degrees of freedom of the 
%approximation basis for each coordinate, r(0) = 1 is not included.
%
%Tiangang Cui, August, 2019

d   = length(ftt.cores);
rs  = ones(1,d);
ns  = ones(1,d);
for k = 1:d-1
    rs(k) = size(ftt.cores{k+1}, 1);
    ns(k) = ftt.oneds{k}.num_nodes;
end
ns(d) = ftt.oneds{d}.num_nodes;

end