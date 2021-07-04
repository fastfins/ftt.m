function [x,x2z] = reference2domain(z, domain)
%Map from the coordinate of the reference domain, z, to the actual nodes, x.
%x2z is the jacobian from x to z.
%
%Tiangang Cui, August, 2019

x2z = 0.5*(domain(:,2)-domain(:,1));
mid = 0.5*(domain(:,2)+domain(:,1));
x   = z(:).*x2z+mid;

x   = reshape(x,size(z));
%x2z = reshape(x2z,size(x,1),[]);

end