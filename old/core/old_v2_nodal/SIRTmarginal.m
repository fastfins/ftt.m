classdef SIRTmarginal < SIRT
    
    methods
        
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
        
        
        function ftt_irt = build_irt_marginal(fttm)
            %Setup data structure used for IRT
            %
            %Inputs:
            %
            %ftt:
            %  A given function tensor train
            %
            %squared:
            %  If the squared mode is turned on, default is yes
            %
            %dir:
            %  the direction of the IRT is constructed,
            %    >0: from the first dimension
            %    <0: from the last dimension
            %
            %Tiangang Cui, August, 2019
            
            if fttm.ng_flag
                
                if fttm.dir == 0
                    error('Not supporting the sampling for marginalisation over centre blocks')
                end
                
                ftt_irt.dir     = fttm.dir;
                ftt_irt.cores   = fttm.cores;
                ftt_irt.oneds   = fttm.oneds;
                ftt_irt.ng_flag = fttm.ng_flag;
                ftt_irt.Y       = fttm.Y;
                
                d = length(fttm.cores);
                % continue building the marginalisation
                ftt_irt.ys = cell(d, 1);
                if fttm.dir > 0
                    % start with the last dim, upper triangular Chol of the mass matrix
                    % ys{d} is built but shouldn't be used
                    Ligeqk  = fttm.Y;
                    for k = d:-1:1
                        nx  = fttm.oneds{k}.num_nodes;
                        rkm = size(fttm.cores{k}, 1);
                        rk  = size(fttm.cores{k}, 3);
                        ftt_irt.ys{k} = reshape( reshape(fttm.cores{k}, rkm*nx, rk)*Ligeqk, rkm, nx, []);
                        B = reshape( oned_mass_r(fttm.oneds{k}, reshape(permute(ftt_irt.ys{k},[2,3,1]),nx,[])), [],rkm);
                        [~, R] = qr(B,0);
                        Ligeqk = R';
                    end
                    ftt_irt.z = sum(sum(Ligeqk.^2, 1));
                else
                    % start with the 1st dim, upper triangular Chol of the mass matrix
                    % ys{1} is built but shouldn't be used
                    Rileqk  = fttm.Y;
                    for k = 1:d
                        nx  = fttm.oneds{k}.num_nodes;
                        rkm = size(fttm.cores{k}, 1);
                        rk  = size(fttm.cores{k}, 3);
                        ftt_irt.ys{k} = reshape( Rileqk*reshape(fttm.cores{k}, rkm,nx*rk), [], nx, rk);
                        B = reshape(oned_mass_r(fttm.oneds{k},reshape(permute(ftt_irt.ys{k},[2,1,3]),nx,[])), [],rk);
                        [~,Rileqk] = qr(B,0);
                    end
                    ftt_irt.z = sum(sum(Rileqk.^2, 1));
                end
                
                ftt_irt.oned_cdfs = cell(size(fttm.oneds));
                for i = 1:length(fttm.oneds)
                    ftt_irt.oned_cdfs{i} = setup_oned_cdf(fttm.oneds{i}, 'squared', true);
                end
                
            else
                % this is just another tensor train
                ftt_irt.dir     = fttm.dir;
                ftt_irt.cores   = fttm.cores;
                ftt_irt.oneds   = fttm.oneds;
                ftt_irt.ng_flag = fttm.ng_flag;
                [ftt_irt.z, ftt_irt.ys]  = build_ftt_cumint(fttm, fttm.dir);
                
                ftt_irt.oned_cdfs = cell(size(fttm.oneds));
                for i = 1:length(fttm.oneds)
                    ftt_irt.oned_cdfs{i} = setup_oned_cdf(fttm.oneds{i}, 'squared', false);
                end
            end
            
        end
        
        function fx = eval_ftt_marginal(data, x)
            %Evaluate the marginalise pdf represented by ftt
            %
            %Inputs:
            %
            %data:
            %  A given marginalised function tensor train
            %
            %x:
            %  input variables
            %
            %Outputs:
            %
            %fx:
            %  marginal density at x
            %
            %Tiangang Cui, August, 2019
            
            if data.ng_flag
                if data.dir > 0
                    % marginalised from the right
                    fxl = eval_ftt_block(data, x, data.dir);
                    fx  = sum((fxl*data.Y).^2, 2)';
                elseif data.dir < 0
                    % marginalised from the left
                    fxg = eval_ftt_block(data, x, data.dir);
                    fx  = sum((data.Y*fxg).^2, 1);
                else
                    % marginalised in the middle
                    nx  = size(x, 2);
                    fxl = eval_ftt_block(data, x(1:data.l_e,:), 1);
                    fxg = eval_ftt_block(data, x(data.r_s:end,:), -1);
                    
                    rkm = size(data.Y, 1);
                    rk  = size(data.Y, 3);
                    tmp = reshape(data.Y, rkm, []);
                    fx  = zeros(1,nx);
                    for i = 1:nx
                        fx(i) = sum((reshape(fxl(i,:)*tmp, [], rk)*fxg(:,i)).^2);
                    end
                end
            else
                fx  = eval_ftt_block(data, x, data.dir);
                if data.dir > 0
                    fx  = fx';
                end
            end
            
            % the normalising constant is not correctly computed
            
        end
        
        function [r,f] = eval_irt_marginal(firt, z)
            %Using IRT to sample from the marginal pdf approximated by firt
            %
            %%%%
            %Inputs:
            %
            %firt:
            %  A data structure built by the function build_irt for evalating the IRT
            %
            %z:
            %  uniform random variables, d x n
            %
            %%%%%
            %Outputs:
            %
            %r:
            %  random variable drawn form the pdf defined by firt
            %
            %f:
            %  pdf function values of those random variables
            %
            %%%%
            %Example:
            %
            %d = 5;
            %a = 0.8;
            %
            %data = setup_ou_process(d, a);
            %func  = @(x) eval_ou_process(data, x); % the pdf
            %poly  = setup_oned(30, 'domain', [-5,5]);
            %
            %opt   = ftt_options('ng_flag', true, 'oned_ref', poly);
            %ftt   = build_ftt(func, d, [], opt);
            %firt  = build_irt(ftt);
            %[r,f] = eval_irt(firt, rand(d, 1E4));
            %
            %
            %Tiangang Cui, August, 2019
            
            %for the non squared version
            if ~firt.ng_flag
                [r,f] = eval_irt(firt, z);
                return;
            end
            
            %for the squared version
            d = length(firt.cores);
            n = size(z,2);
            if size(z,1) ~= d
                error('Input variable must be d x n')
            end
            r = zeros(d,n);
            
            if firt.dir > 0
                frl = ones(n,1);
                for k = 1:d
                    rkm = size(firt.cores{k}, 1);
                    nc  = firt.oned_cdfs{k}.num_nodes;
                    
                    T1  = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
                    pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)';
                    
                    r(k,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(k,:));
                    jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                    ii  = repmat((1:n)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                f   = sum((frl*firt.Y).^2, 2)'/firt.z;
                
            else
                frg = ones(1,n);
                for k = d:-1:1
                    rk  = size(firt.cores{k}, 3);
                    nc  = firt.oned_cdfs{k}.num_nodes;
                    
                    T1  = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
                    pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
                    
                    r(k,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
                    ii  = reshape(1:rk*n, [], 1);
                    jj  = reshape(repmat(1:n, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*n, n);
                    frg = T2'*B;
                end
                f  = sum((firt.Y*frg).^2, 1)/firt.z;
            end
            
        end
        
    end
    
end
