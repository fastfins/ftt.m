function z = int(obj)
% Integrate the entire TT
%

d = length(obj.cores);
if obj.direction > 0
    % last dimenion
    figeqk  = 1;
    for k = d:-1:1
        nx  = obj.oneds{k}.num_nodes;
        rkm = size(obj.cores{k}, 1);
        rk  = size(obj.cores{k}, 3);
        % push the previous dimenion into the current ftt by modifying the
        % coefficients, then ys{k} is used to replace the ftt core at the
        % k-th coordinate for obtaining the integrated function over the
        % last d-k dimensions, ys{k} is a rk-1 by nx matrix
        ys = reshape(reshape(obj.cores{k}, rkm*nx, rk)*figeqk, rkm, nx);
        figeqk  = integral(obj.oneds{k}, ys')';
    end
    z = figeqk;
else
    % first dimenion
    fileqk = 1;
    for k = 1:d
        nx  = obj.oneds{k}.num_nodes;
        rkm = size(obj.cores{k}, 1);
        rk  = size(obj.cores{k}, 3);
        % push the previous dimenion into the current ftt
        % ys{k} is a nx by rk matrix
        ys = reshape(fileqk*reshape(obj.cores{k}, rkm, nx*rk), nx, rk);
        fileqk  = integral(obj.oneds{k}, ys);
    end
    z = fileqk;
end
end