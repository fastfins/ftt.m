function data = build_ftt_marginal(ftt, ind)

d = length(ftt.cores);

if length(ind) ~= length(unique(ind))
    disp('marginalisation indices are not unique')
    ind = unique(ind);
end

if max(ind) > d || min(ind) < 1
    error('Marginalisation indices out of bounds')
end

if length(ind) == d
    error('marginalisation over all indices, should use build_ftt_cumint for integration')
end

data.oneds = ftt.oneds;
data.cores = ftt.cores;
data.dir     = ftt.direction;
data.ng_flag = ftt.ng_flag;
    
if ftt.ng_flag
    ks  = ind(1);
    ke  = ind(end);
    % must marginalise over continous indicies
    if (ke - ks + 1) ~= length(ind)
        error('Can only marginalise over continous indicies in the squared mode')
    end
    
    if ke == d   
        % from the right end
        ind = fliplr(ind);    
        R   = 1;
        for i = 1:length(ind)
            k   = ind(i);
            nx  = data.oneds{k}.num_nodes;
            rkm = size(data.cores{k}, 1);
            rk  = size(data.cores{k}, 3);
            
            tmp = reshape( reshape(data.cores{k}, rkm*nx, rk)*R', rkm, nx, []);
            B   = reshape( oned_mass_r(data.oneds{k}, reshape(permute(tmp,[2,3,1]),nx,[])), [],rkm);
            [~,R] = qr(B,0);
        end
        data.Y   = R';
        % marginalised from the right (>0) 
        data.dir = 1;
    elseif ks == 1
        R   = 1;
        for i = 1:length(ind)
            k   = ind(i);
            nx  = data.oneds{k}.num_nodes;
            rkm = size(data.cores{k}, 1);
            rk  = size(data.cores{k}, 3);
            
            tmp = reshape( R*reshape(data.cores{k}, rkm,nx*rk), [], nx, rk);
            B   = reshape(oned_mass_r(data.oneds{k},reshape(permute(tmp,[2,1,3]),nx,[])), [],rk);
            [~,R] = qr(B,0);
        end
        data.Y   = R;
        % marginalised from the left (<0) 
        data.dir = -1;
    else
        disp('Marginalisation over centre blocks is very costly to operate, try to avoid')
        % in the middle
        rl2 = size(ftt.cores{ks-1}, 3);
        rr1 = size(ftt.cores{ke+1}, 1);
        Rc  = 1;
        for i = 1:length(ind)
            k   = ind(i);
            nx  = data.oneds{k}.num_nodes;
            rkm = size(data.cores{k}, 1);
            rk  = size(data.cores{k}, 3);
            
            R   = oned_mass_r(data.oneds{k}, reshape(permute(data.cores{k}, [2,1,3]), nx,rkm*rk));
            
            tmp = Rc*reshape(permute(reshape(R, nx,rkm,rk),[2,1,3]), rkm, []);
            tmp = reshape(permute(reshape(tmp, [], rl2, nx,rk), [1,3,2,4]), [], rl2*rk);
            [~,Rc]  = qr(tmp, 0);
            Rc  = reshape(Rc, [], rk);
        end
        
        data.Y  = permute(reshape(Rc, [], rl2, rr1), [2,1,3]);
        
        % marginalised in the middle (=0) 
        data.dir = 0;
        data.l_e = ks-1;
        data.r_s = ks;
    end
    
    % delete marginalised blocks
    data.oneds(ind) = [];
    data.cores(ind) = [];
    
    % rebuilt the squared root ftt
    
else
    %need to identify continuous block from the right end
    [ind1, ind2] = process_indices(d, ind);
    
    %the left blocks
    if ~isempty(ind1)
        for i = 1:length(ind1)
            k   = ind1(i);
            nx  = data.oneds{k}.num_nodes;
            
            rkm = size(data.cores{k}, 1);
            rk  = size(data.cores{k}, 3);
            tmp = oned_integral(data.oneds{k}, reshape(permute(data.cores{k}, [2,1,3]), nx,rkm*rk));
            tmp = reshape(tmp, rkm, rk);
            data.cores{k+1} = reshape(tmp*reshape(data.cores{k+1}, rk, []), rkm, nx, []);
        end
    end
    
    % the right blocks
    if ~isempty(ind2)
        for i = 1:length(ind2)
            k = ind2(i);
            nx  = data.oneds{k}.num_nodes;
            rkm = size(data.cores{k}, 1);
            rk  = size(data.cores{k}, 3); % rk should always be 1
            tmp = oned_integral(data.oneds{k}, reshape(permute(data.cores{k}, [2,1,3]), nx,rkm*rk));
            tmp = reshape(tmp, rkm, rk);
            data.cores{k-1} = reshape(reshape(data.cores{k-1}, [], rkm)*tmp, [], nx, rk);
        end
    end
    
    % delete marginalised blocks
    data.oneds(ind) = [];
    data.cores(ind) = [];
end

    
end

%{
% in the middle
M   = 1;
rl2 = size(ftt.cores{ks-1}, 3);
rr1 = size(ftt.cores{ke+1}, 1);
Rc  = 1;
for i = 1:length(ind)
    k   = ind(i);
    nx  = data.oneds{k}.num_nodes;
    rkm = size(data.cores{k}, 1);
    rk  = size(data.cores{k}, 3);

    R   = oned_mass_r(data.oneds{k}, reshape(permute(data.cores{k}, [2,1,3]), nx,rkm*rk));
    Mc  = reshape(permute(reshape(R'*R, rkm, rk, rkm, rk), [1 3 2 4]), rkm^2, rk^2);
    M   = M*Mc;

    tmp = Rc*reshape(permute(reshape(R, nx,rkm,rk),[2,1,3]), rkm, []);
    tmp = reshape(permute(reshape(tmp, [], rl2, nx,rk), [1,3,2,4]), [], rl2*rk);
    [~,Rc] = qr(tmp, 0);
    Rc  = reshape(Rc, [], rk);
end
%nxl = ftt.oneds{ks-1}.num_nodes;
%rl1 = size(ftt.cores{ks-1}, 1);
rl2 = size(ftt.cores{ks-1}, 3);
%nxr = ftt.oneds{ke+1}.num_nodes;
rr1 = size(ftt.cores{ke+1}, 1);
%rr2 = size(ftt.cores{ke+1}, 3);

% M should be algebraically symmetric and U, V should be algebraically same
% M   = reshape(permute(reshape(M,rl2,rl2,rr1,rr1), [1 3 2 4]), rl2*rr1, rl2*rr1);
%{
[L,D] = ldl(0.5*(M+M'));
d  = diag(D);
ii = d > 0;
data.Y   = permute(reshape(L(:,ii)*diag(d(ii).^(0.5)), rl2, rr1, []), [1,3,2]);
%}

%{
[U,D]   = eig(0.5*(M+M'));
[d,jj]  = sort(diag(D), 'descend');
s   = d.^(0.5);
ii  = (d > eps);
data.Y   = permute(reshape(U(:,jj(ii))*diag(s(ii)), rl2, rr1, []), [1,3,2]);
%}        

Rc  = reshape(Rc, [], rl2*rr1);


% marginalised in the middle (=0) 
data.dir = 0;
data.l_e = ks-1;
data.r_s = ks;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind1, ind2] = process_indices(d, ind)

ind = reshape(ind, 1, []);

if ind(end) == d
    jr = d;
    for i = (length(ind)-1):-1:1
        j = ind(i);
        if j+1 == jr
            % if the index is continous, move to the left
            jr = j;
        else
            % otherwise stop searching, jr is the last stopping index
            % jr:d are continous block
            break;
        end
    end
    ind2 = jr:d;
    ind1 = setdiff(ind, ind2);
    ind2 = fliplr(ind2);
else
    ind1 = ind;
    ind2 = [];
end

end