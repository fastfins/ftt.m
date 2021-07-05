function obj = round(obj, thres)
% round the TT cores
%
if nargin == 1
    thres = obj.opt.local_tol;
end
[d,rs,~] = size(obj);
%
obj.direction = -obj.direction;
if obj.direction > 0
    ind = 1:(d-1);
else
    ind = d:-1:2;
end
%
% start
for k = ind
    if obj.direction > 0
        if k == 1
            Jx_left = [];
        else
            Jx_left = obj.interp_x{k-1};
        end
        [obj.cores{k}, obj.interp_x{k}, obj.cores{k+1}] = FTT.build_basis_svd(obj.oneds{k},...
            Jx_left, obj.cores{k+1}, obj.cores{k}, ...
            obj.direction, obj.opt.int_method, thres, obj.opt.max_rank);
        rs(k) = size(obj.cores{k}, 3);
    else
        if k == d
            Jx_right = [];
        else
            Jx_right = obj.interp_x{k+1};
        end
        [obj.cores{k}, obj.interp_x{k}, obj.cores{k-1}] = FTT.build_basis_svd(obj.oneds{k},...
            Jx_right, obj.cores{k-1}, obj.cores{k}, ...
            obj.direction, obj.opt.int_method, thres, obj.opt.max_rank);
        rs(k-1) = size(obj.cores{k}, 1);
    end
end

if strcmp (obj.opt.tt_method, 'amen')
    obj.res_w = [];
end

disp('rounded TT')
disp(rs)
end