function [z,x2z] = domain2reference(x, domain)
%Map from the actual nodes, x, to the coordinate of the reference domain, z. 
%x2z is the jacobian from x to z.
%
%Tiangang Cui, August, 2019

x2z = 0.5*(domain(:,2)-domain(:,1));
mid = 0.5*(domain(:,2)+domain(:,1));
z   = (x(:)-mid)./x2z;

z   = reshape(z,size(x));
%x2z = reshape(x2z,size(x,1),[]);

end