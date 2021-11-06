function [gx,fx] = grad(obj, x)
% Evaluate the gradient of the TT
%   [g,f] = GRAD(tt, x)
%
%   x   - input variables, k x n.
%   g   - gradient at x, k x n
%   f   - function values at x, 1 x n

nx  = size(x, 2);
fxl = cell(obj.d,1);
fxr = cell(obj.d,1);
    
for j = 1:obj.d
    nj  = size(obj.cores{j}, 2);
    rjm = size(obj.cores{j}, 1);
    %
    tmp = reshape(permute(obj.cores{j}, [2,1,3]), nj, []);
    % rjm nx rj
    if j == 1
        fxl{j} = reshape(permute(reshape( eval(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
    else
        T   = reshape(permute(reshape( eval(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
        % how to speed up this part?
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        Bf  = sparse(ii(:), jj(:), fxl{j-1}(:), nx, rjm*nx);
        %
        fxl{j}  = Bf*T; % nx by rj
    end
end
%
for j = obj.d:-1:1
    nj  = size(obj.cores{j}, 2);
    rj  = size(obj.cores{j}, 3);
    %
    tmp = reshape(permute(obj.cores{j}, [2,3,1]), nj, []);
    % rj nx rjm
    if j == obj.d
        fxr{j} = reshape(permute(reshape( eval(obj.oneds{j},tmp,x(j,:)), nx, rj, []), [3,1,2]), rjm, []);
    else
        T = reshape(permute(reshape( eval(obj.oneds{j},tmp,x(j,:)), nx, rj, []), [2,1,3]), rj*nx, []);
        % how to speed up this part?
        ii  = reshape(1:rj*nx, [], 1);
        jj  = reshape(repmat(1:nx, rj, 1), [], 1);
        Bf  = sparse(ii, jj, fxr{j+1}(:), rj*nx, nx);
        %
        fxr{j} = T'*Bf; % rjm by nx
    end
end
%
gx = zeros(obj.d,nx);
for j = 1:obj.d
    nj  = size(obj.cores{j}, 2);
    rjm = size(obj.cores{j}, 1);
    %
    tmp = reshape(permute(obj.cores{j}, [2,1,3]), nj, []);
    % (rjm nx) by rj
    D = reshape(permute(reshape( eval_deri(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
    if j == 1
        gx(j,:) = sum(D'.*fxr{j+1}, 1);
    elseif j == obj.d
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        % nx by (rjm nx)
        Bl  = sparse(ii(:), jj(:), fxl{j-1}(:), nx, rjm*nx);
        % fxr{j+1}: rj by nx,   Bl*D: nx by rj
        gx(j,:) =(Bl*D)';
    else
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        % nx by (rjm nx)
        Bl  = sparse(ii(:), jj(:), fxl{j-1}(:), nx, rjm*nx);
        % fxr{j+1}: rj by nx,   Bl*D: nx by rj
        gx(j,:) = sum((Bl*D)'.*fxr{j+1}, 1);
    end
end
fx = fxr{1};

if obj.opt.sqrt_flag
    fx = fx.^2;
    gx = 2*(gx.*fx);
end

end

%{
function [gx,fx] = grad(obj, x)
% Evaluate the gradient of the TT
%   [g,f] = GRAD(tt, x)
%
%   x   - input variables, k x n.
%   g   - gradient at x, k x n
%   f   - function values at x, 1 x n

nx = size(x, 2);
Ts = cell(obj.d,1);
Ds = cell(obj.d,1);

if obj.direction > 0
    for j = 1:obj.d
        nj  = size(obj.cores{j}, 2);
        rjm = size(obj.cores{j}, 1);
        %
        if j < obj.d || (size(obj.cores{j}, 3) > 1 && size(obj.cores{j}, 4) == 1)
            tmp = reshape(permute(obj.cores{j}, [2,1,3]), nj, []);
        else
            tmp = reshape(permute(reshape(obj.cores{j}, rjm, nj, []), [2,1,3]), nj, []);
        end
        % rjm nx rj
        Ts{j} = permute(reshape( eval(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]);
        Ds{j} = permute(reshape( eval_deri(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]);
    end
    %
    fx = ones(nx,1);
    gx = ones(nx,obj.d);
    for j = 1:obj.d
        rjm = size(obj.cores{j}, 1);
        rj  = size(obj.cores{j}, 3);
        m   = size(obj.cores{j}, 4);
        %
        T   = reshape(Ts{j}, rjm*nx, []);
        D   = reshape(Ds{j}, rjm*nx, []);
        % how to speed up this part?
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        Bf  = sparse(ii(:), jj(:), fx(:), nx, rjm*nx);
        %
        fx  = Bf*T; % nx by rj*m
        %
        tgx = zeros(nx*rj*m,obj.d);
        for k = 1:obj.d
            Bg  = sparse(ii(:), jj(:), gx(:,k), nx, rjm*nx);
            if k == j
                tgx(:,k) = reshape(Bg*D, [], 1);
            else
                tgx(:,k) = reshape(Bg*T, [], 1);
            end
        end
        gx = tgx;
    end
    fx = fx';
    gx = gx';
else
    for j = obj.d:-1:1
        nj  = size(obj.cores{j}, 2);
        rj  = size(obj.cores{j}, 3);
        %
        if j > 1 || (size(obj.cores{j}, 1) > 1 && size(obj.cores{j}, 4) == 1)
            tmp = reshape(permute(obj.cores{j}, [2,3,1]), nj, []);
        else
            % collapse the first dimension = 1
            tmp = reshape(obj.cores{j}, nj, []);
        end
        % rj nx rjm
        Ts{j} = permute(reshape( eval(obj.oneds{j},tmp,x(j,:)), nx, rj, []), [2,1,3]);
        Ds{j} = permute(reshape( eval_deri(obj.oneds{j},tmp,x(j,:)), nx, rj, []), [2,1,3]);
    end
    %
    fx = ones(1,nx);
    gx = ones(obj.d,nx);
    for j = obj.d:-1:1
        rjm = size(obj.cores{j}, 1);
        rj  = size(obj.cores{j}, 3);
        m   = size(obj.cores{j}, 4);
        %
        T   = reshape(Ts{j}, rj*nx, []);
        D   = reshape(Ds{j}, rj*nx, []);
        % how to speed up this part?
        ii  = reshape(1:rj*nx, [], 1);
        jj  = reshape(repmat(1:nx, rj, 1), [], 1);
        Bf  = sparse(ii, jj, fx(:), rj*nx, nx);
        %
        fx  = T'*Bf; % rjm by nx
        %
        tgx = zeros(obj.d,rjm*nx*m);
        for k = 1:obj.d
            Bg  = sparse(ii(:), jj(:), gx(k,:)', rj*nx, nx);
            if k == j
                tgx(k,:) = reshape(D'*Bg, [], 1)';
            else
                tgx(k,:) = reshape(T'*Bg, [], 1)';
            end
        end
        gx = tgx;
    end
end

if obj.opt.sqrt_flag
    fx = fx.^2;
    gx = 2*(gx.*fx);
end

end

%}
