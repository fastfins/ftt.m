classdef SIRT
    
    properties
        dir
        cores
        oneds
        sqrt_flag
    end
    
    methods (Static)
        function T = eval_oned_core_213(oned, core, x)
            
            rkm = size(core, 1);
            nn  = oned.num_nodes;
            nx  = length(x);
            
            % evaluate the updated basis function
            tmp = eval_oned_nodes2x(oned, x(:), reshape(permute(core, [2,1,3]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
            
        end
        
        function T = eval_oned_core_213_deri(oned, core, x)
            
            rkm = size(core, 1);
            nn  = oned.num_nodes;
            nx  = length(x);
            
            % evaluate the updated basis function
            tmp = eval_oned_nodes2x_deri(oned, x(:), reshape(permute(core, [2,1,3]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
            
        end
        
        function T = eval_oned_core_231(oned, core, x)
            
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            
            % evaluate the updated basis function
            tmp = eval_oned_nodes2x(oned, x(:), reshape(permute(core, [2,3,1]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
            
        end
        
        function T = eval_oned_core_231_deri(oned, core, x)
            
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            
            % evaluate the updated basis function
            tmp = eval_oned_nodes2x_deri(oned, x(:), reshape(permute(core, [2,3,1]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
            
        end
        
        
    end
    
    methods
        function obj = IRT(ftt, varargin)
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
            
            defaultDir      = 1;
            
            p = inputParser;
            %
            addRequired (p,'ftt');
            addParameter(p,'dir',    defaultDir,    @(x) isnumeric(x) && isscalar(x) && (x~=0));
            
            parse(p,ftt,varargin{:});
            
            dir     = p.Results.dir;
            
            obj.dir     = dir;
            obj.cores   = ftt.cores;
            obj.oneds   = ftt.oneds;
            obj.sqrt_flag = ftt.sqrt_flag;
            [ftt_irt.z, ftt_irt.ys]  = cumint(ftt, dir);
            
            ftt_irt.oned_cdfs = cell(size(ftt.oneds));
            for i = 1:length(ftt.oneds)
                ftt_irt.oned_cdfs{i} = setup_oned_cdf(ftt.oneds{i}, 'squared', ftt.ng_flag);
            end
        end
        
        
        
        function [z,ys] = build_ftt_cumint(ftt, dir)
            %Marginalise the pdf represented by ftt sequentially dimension by dimension
            %
            %Inputs:
            %
            %ftt:
            %  A given function tensor train
            %
            %dir:
            %  the direction of the marginalisation
            %    >0: from the last dimension, marginalise to the first
            %    <0: from the first dimension, marginalise to the last
            %
            %Outputs:
            %
            %z:
            %  Normalising constant
            %
            %ys:
            %  A cell array holding the marginalised coefficents for marginalisation
            %  up to each dimension
            %
            %Tiangang Cui, August, 2019
            
            d  = length(ftt.cores);
            ys = cell(d, 1);
            
            if ftt.ng_flag
                
                if dir > 0
                    % start with the last dim, upper triangular Chol of the mass matrix
                    % ys{d} is built but shouldn't be used
                    Ligeqk  = 1;
                    for k = d:-1:1
                        nx  = ftt.oneds{k}.num_nodes;
                        rkm = size(ftt.cores{k}, 1);
                        rk  = size(ftt.cores{k}, 3);
                        % push the previous dimenion into the current ftt by modifying the
                        % coefficients, then ys{k} is used to replace the ftt core at the
                        % k-th coordinate for obtaining the integrated function over the
                        % last d-k dimensions, ys{k} is a rk-1 by nx by rk tensor, as we
                        % still need rk functions to evaluate the function
                        ys{k} = reshape( reshape(ftt.cores{k}, rkm*nx, rk)*Ligeqk, rkm, nx, rk);
                        
                        %disp(size(ftt.cores{k}))
                        %disp(size(ys{k}))
                        
                        B = reshape( oned_mass_r(ftt.oneds{k}, reshape(permute(ys{k},[2,3,1]),nx,[])), [],rkm);
                        [~, R] = qr(B,0);
                        if size(R, 1) < size(R, 2)
                            R = R(:,1:size(R, 1));
                        end
                        Ligeqk = R';
                    end
                    z = sum(sum(Ligeqk.^2, 1));
                    %Ligeqk
                else
                    % start with the 1st dim, upper triangular Chol of the mass matrix
                    % ys{1} is built but shouldn't be used
                    Rileqk  = 1;
                    for k = 1:d
                        nx  = ftt.oneds{k}.num_nodes;
                        rkm = size(ftt.cores{k}, 1);
                        rk  = size(ftt.cores{k}, 3);
                        % push the previous dimenion into the current ftt by modifying the
                        % coefficients, then ys{k} is used to replace the ftt core at the
                        % k-th coordinate for obtaining the integrated function over the
                        % last d-k dimensions, ys{k} is a rk-1 by nx by rk tensor, as we
                        % still need rk-1 functions to evaluate the function
                        ys{k} = reshape( Rileqk*reshape(ftt.cores{k}, rkm,nx*rk), [], nx, rk);
                        
                        %disp(size(ftt.cores{k}))
                        %disp(size(ys{k}))
                        
                        B = reshape(oned_mass_r(ftt.oneds{k},reshape(permute(ys{k},[2,1,3]),nx,[])), [],rk);
                        [~,Rileqk] = qr(B,0);
                        %if size(Rileqk, 1) < size(Rileqk, 2)
                        %   Rileqk = Rileqk(:,1:size(Rileqk, 1));
                        %end
                        
                        %size(B)
                        %size(Rileqk)
                        
                    end
                    z = sum(sum(Rileqk.^2, 1));
                    %Rileqk
                end
            else
                if dir > 0
                    % last dimenion
                    figeqk  = 1;
                    for k = d:-1:1
                        nx  = ftt.oneds{k}.num_nodes;
                        rkm = size(ftt.cores{k}, 1);
                        rk  = size(ftt.cores{k}, 3);
                        % push the previous dimenion into the current ftt by modifying the
                        % coefficients, then ys{k} is used to replace the ftt core at the
                        % k-th coordinate for obtaining the integrated function over the
                        % last d-k dimensions, ys{k} is a rk-1 by nx matrix
                        ys{k}   = reshape(reshape(ftt.cores{k}, rkm*nx, rk)*figeqk, rkm, nx);
                        figeqk  = oned_integral(ftt.oneds{k}, ys{k}')';
                    end
                    z   = figeqk;
                else
                    % first dimenion
                    fileqk = 1;
                    for k = 1:d
                        nx  = ftt.oneds{k}.num_nodes;
                        rkm = size(ftt.cores{k}, 1);
                        rk  = size(ftt.cores{k}, 3);
                        % push the previous dimenion into the current ftt
                        % ys{k} is a nx by rk matrix
                        ys{k} = reshape(fileqk*reshape(ftt.cores{k}, rkm, nx*rk), nx, rk);
                        fileqk  = oned_integral(ftt.oneds{k}, ys{k});
                    end
                    z   = fileqk;
                    % end
                end
                
            end
            
        end
        
        function [r,f] = eval_irt(firt, z)
            %Using IRT to sample from the target pdf approximated by firt
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
                    
                    if firt.ng_flag
                        T1  = reshape( eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
                        pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)';
                    else
                        pk  = firt.ys{k}'*frl';
                    end
                    
                    r(k,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(k,:));
                    jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                    ii  = repmat((1:n)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                
                if firt.ng_flag
                    f   = frl'.^2/firt.z;
                else
                    f   = frl'/firt.z;
                end
                
            else
                frg = ones(1,n);
                for k = d:-1:1
                    rk  = size(firt.cores{k}, 3);
                    nc  = firt.oned_cdfs{k}.num_nodes;
                    
                    if firt.ng_flag
                        T1  = reshape( eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
                        pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
                    else
                        pk  = firt.ys{k}*frg;
                    end
                    
                    r(k,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
                    ii  = reshape(1:rk*n, [], 1);
                    jj  = reshape(repmat(1:n, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*n, n);
                    frg = T2'*B;
                end
                
                if firt.ng_flag
                    f   = frg.^2/firt.z;
                else
                    f   = frg/firt.z;
                end
                
            end
            
        end
        
        function [r,f] = eval_cond_irt(firt, x, z)
            %Using IRT to draw conditional samples from the target pdf approximated by
            %firt
            %
            %%%%
            %Inputs:
            %
            %firt:
            %  A data structure built by the function build_irt for evalating the IRT
            %
            %z:
            %  uniform random variables, m x n, used for generating conditional samples
            %
            %x:
            %  given random variable: (d-m) x n
            %
            %%%%%
            %Outputs:
            %
            %r:
            %  random variable drawn form the pdf p(.|x) defined by firt
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
            
            
            d  = length(firt.cores);
            nr = size(z,2);
            dr = size(z,1);
            nx = size(x,2);
            dx = size(x,1);
            
            if dx == 0 || dr == 0
                error('dimension of x or the dimension of z should be nonzero')
            end
            
            if dr + dx ~= d
                error('dimension of x and the dimension of z mismatch the dimension of the joint')
            end
            
            if nx ~= nr
                if nx ~= 1
                    error('number of x and the number of z mismatch')
                end
            end
            
            r = zeros(size(z));
            
            % first compute the marginal
            if firt.dir > 0
                %from the last dimension, marginalise to the first
                %order of the samples (x, r)
                %first evaluate the marginal for x
                
                frl = ones(nr,1);
                for k = 1:(dx-1)
                    rkm = size(firt.cores{k}, 1);
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, x(k,:));
                    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
                    ii  = repmat((1:nx)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                
                k   = dx;
                rkm = size(firt.cores{k}, 1);
                % evaluate the updated basis function
                % for the marginal
                if firt.ng_flag
                    T2  = eval_oned_core_213(firt.oneds{k}, firt.ys{k}, x(k,:));
                    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
                    ii  = repmat((1:nx)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
                    frl_m   = B*T2;
                    fm  = sum(frl_m.^2, 2)';
                else
                    fxk = eval_oned_nodes2x(firt.oneds{k}, x(k,:), firt.ys{k}');
                    % mutiply together
                    fm  = sum(frl.*fxk, 2)';
                end
                
                %for the frl, evaluate the updated basis function
                T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, x(k,:));
                jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
                ii  = repmat((1:nx)', 1, rkm);
                B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
                frl = B*T2;
                
                %then generate the conditional samples
                for j = 1:dr
                    k = dx+j;
                    rkm = size(firt.cores{k}, 1);
                    nc  = firt.oned_cdfs{k}.num_nodes;
                    
                    %bk  = compute_basis(firt.oneds{k}, firt.oned_cdfs{k});
                    if firt.ng_flag
                        T1  = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
                        pk  = reshape(sum(reshape(frl*T1, nr*nc, []).^2, 2), nr, nc)';
                    else
                        pk  = firt.ys{k}'*frl';
                    end
                    
                    r(j,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(j,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(j,:));
                    jj  = reshape(reshape(1:rkm*nr, rkm, nr)', [], 1);
                    ii  = repmat((1:nr)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), nr, rkm*nr);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                
                if firt.ng_flag
                    f   = frl'.^2;
                else
                    f   = frl';
                end
            else
                %from the first dimension, marginalise to the last
                %order of the samples (r, x)
                %first evaluate the marginal for x
                
                frg = ones(1,nr);
                for j = dx:-1:2
                    k   = dr + j;
                    rk  = size(firt.cores{k}, 3);
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, x(j,:));
                    ii  = reshape(1:rk*nx, [], 1);
                    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*nx, nx);
                    frg = T2'*B;
                end
                
                k   = dr+1;
                rk  = size(firt.cores{k}, 3);
                % evaluate the updated basis function
                
                % for marginal
                if firt.ng_flag
                    T2  = eval_oned_core_231(firt.oneds{k}, firt.ys{k}, x(1,:));
                    ii  = reshape(1:rk*nx, [], 1);
                    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*nx, nx);
                    frg_m   = T2'*B;
                    fm  = sum(frg_m.^2, 1);
                else
                    fxk = eval_oned_nodes2x(firt.oneds{k}, x(1,:), firt.ys{k});
                    % mutiply together
                    fm  = sum(fxk'.*frg, 1);
                end
                
                % for frg
                T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, x(1,:));
                ii  = reshape(1:rk*nx, [], 1);
                jj  = reshape(repmat(1:nx, rk, 1), [], 1);
                B   = sparse(ii, jj, frg(:), rk*nx, nx);
                frg = T2'*B;
                
                %then generate the conditional samples
                for k = dr:-1:1
                    rk  = size(firt.cores{k}, 3);
                    nc  = firt.oned_cdfs{k}.num_nodes;
                    
                    if firt.ng_flag
                        T1  = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
                        pk  = reshape(sum(reshape(T1*frg, [], nc*nr).^2, 1), nc, nr);
                    else
                        pk  = firt.ys{k}*frg;
                    end
                    
                    r(k,:)  = sample_oned_cdf(firt.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
                    ii  = reshape(1:rk*nr, [], 1);
                    jj  = reshape(repmat(1:nr, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*nr, nr);
                    frg = T2'*B;
                end
                
                if firt.ng_flag
                    f   = frg.^2;
                else
                    f   = frg;
                end
            end
            
            f = f./fm;
            
        end
        
        function z = eval_rt(firt, r)
            %Using IRT to sample from the target pdf approximated by firt
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
            
            
            d = length(firt.cores);
            n = size(r,2);
            if size(r,1) ~= d
                error('Input variable must be d x n')
            end
            z = zeros(d,n);
            
            if firt.dir > 0
                frl = ones(n,1);
                for k = 1:d
                    rkm = size(firt.cores{k}, 1);
                    nc  = firt.oned_cdfs{k}.num_nodes;
                    
                    if firt.ng_flag
                        T1  = reshape( eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
                        pk  = reshape(sum(reshape(frl*T1, n*nc, []).^2, 2), n, nc)';
                    else
                        pk  = firt.ys{k}'*frl';
                    end
                    
                    z(k,:)  = eval_oned_cdf(firt.oned_cdfs{k}, pk, r(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(k,:));
                    jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                    ii  = repmat((1:n)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                
            else
                frg = ones(1,n);
                for k = d:-1:1
                    rk  = size(firt.cores{k}, 3);
                    nc  = firt.oned_cdfs{k}.num_nodes;
                    
                    if firt.ng_flag
                        T1  = reshape( eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
                        pk  = reshape(sum(reshape(T1*frg, [], nc*n).^2, 1), nc, n);
                    else
                        pk  = firt.ys{k}*frg;
                    end
                    
                    z(k,:)  = eval_oned_cdf(firt.oned_cdfs{k}, pk, r(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
                    ii  = reshape(1:rk*n, [], 1);
                    jj  = reshape(repmat(1:n, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*n, n);
                    frg = T2'*B;
                end
                
            end
            
        end
        
        function J = eval_rt_jac(firt, r, z)
            %Using IRT to sample from the target pdf approximated by firt
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
            
            
            d = length(firt.cores);
            n = size(r,2);
            if size(r,1) ~= d
                error('Input variable must be d x n')
            end
            J = zeros(d,n*d);
            
            if firt.dir > 0
                % function value of each block
                block_ftt = cell(d,1);
                block_mar = cell(d,1);
                T = cell(d,1);
                block_ftt_d = cell(d,1);
                for k = 1:d
                    block_ftt{k} = eval_oned_core_213(firt.oneds{k}, firt.cores{k}, r(k,:));
                    block_mar{k} = eval_oned_core_213(firt.oneds{k}, firt.ys{k}, r(k,:));
                    %
                    rkm  = size(firt.cores{k}, 1);
                    T{k} = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), rkm, []);
                    %
                    block_ftt_d{k} = eval_oned_core_213_deri(firt.oneds{k}, firt.cores{k}, r(k,:));
                end
                %
                F = cell(d,1); % accumulated ftt
                G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
                F{1} = block_ftt{1};
                G{1} = block_mar{1};
                for k = 2:d
                    rkm = size(firt.cores{k}, 1);
                    jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                    ii  = repmat((1:n)', 1, rkm);
                    B   = sparse(ii(:), jj(:), F{k-1}(:), n, rkm*n);
                    F{k} = B*block_ftt{k};
                    G{k} = B*block_mar{k};
                end
                %
                Fm = cell(d,1); % accumulated ftt
                for k = 1:d
                    Fm{k} = sum(G{k}.^2, 2)';
                end
                %
                % j is the index of the coordinate of differentiation
                for j = 1:d
                    ind = (1:d:n*d) + (j-1);
                    %
                    % the diagonal
                    if j == 1
                        J(1,ind) = Fm{1}/firt.z;
                    else
                        J(j,ind) = Fm{j}./Fm{j-1};
                    end
                    %
                    % derivative of the FTT
                    drl = block_ftt_d{j};
                    % derivative of the FTT, 2nd term, for the d(j+1)/dj term
                    mrl = eval_oned_core_213_deri(firt.oneds{j}, firt.ys{j}, r(j,:));
                    if j > 1
                        rjm = size(firt.cores{j}, 1);
                        jj  = reshape(reshape(1:rjm*n, rjm, n)', [], 1);
                        ii  = repmat((1:n)', 1, rjm);
                        B   = sparse(ii(:), jj(:), F{j-1}(:), n, rjm*n);
                        drl = B*drl;
                        mrl = B*mrl;
                    end
                    if j < d
                        J(j+1,ind) = J(j+1,ind) - sum(G{j}.*mrl,2)'.*z(j+1,:);
                    end
                    %
                    for k = (j+1):d
                        % accumulate the j-th block and ealuate the integral
                        rkm = size(firt.cores{k}, 1);
                        nc  = firt.oned_cdfs{k}.num_nodes;
                        pk  = reshape(sum( reshape((F{k-1}*T{k}).*(drl*T{k}), n*nc, []), 2), n, nc)';
                        % the first term
                        J(k,ind) = J(k,ind) + reshape(eval_oned_cdf_deri(firt.oned_cdfs{k}, pk, r(k,:)), 1, []);
                        %
                        if k < d
                            jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                            ii  = repmat((1:n)', 1, rkm);
                            B   = sparse(ii(:), jj(:), drl(:), n, rkm*n);
                            % the second term, for the d(k+1)/dj term
                            mrl = B*block_mar{k};
                            J(k+1,ind) = J(k+1,ind) - sum(G{k}.*mrl,2)'.*z(k+1,:);
                            % acumulate
                            drl = B*block_ftt{k};
                        end
                        J(k,ind) = 2*J(k,ind)./Fm{k-1};
                    end
                    %
                end
                
            else
                % function value of each block
                block_ftt = cell(d,1);
                block_mar = cell(d,1);
                T = cell(d,1);
                block_ftt_d = cell(d,1);
                for k = 1:d
                    block_ftt{k} = eval_oned_core_231(firt.oneds{k}, firt.cores{k}, r(k,:));
                    block_mar{k} = eval_oned_core_231(firt.oneds{k}, firt.ys{k}, r(k,:));
                    rk   = size(firt.cores{k}, 3);
                    T{k} = reshape(eval_oned_core_213(firt.oneds{k}, firt.ys{k}, firt.oned_cdfs{k}.nodes(:)), [], rk);
                    block_ftt_d{k} = eval_oned_core_231_deri(firt.oneds{k}, firt.cores{k}, r(k,:));
                end
                %
                F = cell(d,1); % accumulated ftt
                G = cell(d,1); % sum( G^.2 ) is the marginal, each G{k} is nxr
                F{d} = block_ftt{d}';
                G{d} = block_mar{d}';
                for k = d-1:-1:1
                    rk  = size(firt.cores{k}, 3);
                    ii  = reshape(1:rk*n, [], 1);
                    jj  = reshape(repmat(1:n, rk, 1), [], 1);
                    B   = sparse(ii, jj, F{k+1}(:), rk*n, n);
                    F{k} = block_ftt{k}'*B;
                    G{k} = block_mar{k}'*B;
                end
                %
                Fm = cell(d,1); % accumulated ftt
                for k = 1:d
                    Fm{k} = sum(G{k}.^2, 1);
                end
                %
                % j is the index of the coordinate of differentiation
                for j = d:-1:1
                    ind = (1:d:n*d) + (j-1);
                    %
                    % the diagonal
                    if j == d
                        J(d,ind) = Fm{d}/firt.z;
                    else
                        J(j,ind) = Fm{j}./Fm{j+1};
                    end
                    %
                    % derivative of the FTT
                    drg = block_ftt_d{j}';
                    % derivative of the FTT, 2nd term, for the d(j+1)/dj term
                    mrg = eval_oned_core_231_deri(firt.oneds{j}, firt.ys{j}, r(j,:))';
                    if j < d
                        rj  = size(firt.cores{j}, 3);
                        ii  = reshape(1:rj*n, [], 1);
                        jj  = reshape(repmat(1:n, rj, 1), [], 1);
                        B   = sparse(ii, jj, F{j+1}(:), rj*n, n);
                        drg = drg*B;
                        mrg = mrg*B;
                    end
                    if j > 1
                        J(j-1,ind) = J(j-1,ind) - sum(G{j}.*mrg,1).*z(j-1,:);
                    end
                    %
                    for k = (j-1):-1:1
                        % accumulate the j-th block and ealuate the integral
                        rk  = size(firt.cores{k}, 3);
                        nc  = firt.oned_cdfs{k}.num_nodes;
                        pk  = reshape(sum(reshape((T{k}*F{k+1}).*(T{k}*drg), [], nc*n), 1), nc, n);
                        % the first term
                        J(k,ind) = J(k,ind) + reshape(eval_oned_cdf_deri(firt.oned_cdfs{k}, pk, r(k,:)), 1, []);
                        %
                        if k > 1
                            ii  = reshape(1:rk*n, [], 1);
                            jj  = reshape(repmat(1:n, rk, 1), [], 1);
                            B   = sparse(ii, jj, drg(:), rk*n, n);
                            % the second term, for the d(k+1)/dj term
                            mrg = block_mar{k}'*B;
                            J(k-1,ind) = J(k-1,ind) - sum(G{k}.*mrg,1).*z(k-1,:);
                            % acumulate
                            drg = block_ftt{k}'*B;
                        end
                        J(k,ind) = 2*J(k,ind)./Fm{k+1};
                    end
                    %
                end
            end
            
        end
        
    end
    
end