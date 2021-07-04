classdef IRT < IRTcommon
    
    properties
        dir
        ftt
        ys
        z
        oned_cdfs
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
            
            defaultDir      = 1;
            p = inputParser;
            %
            addRequired (p,'ftt');
            addOptional(p,'dir', defaultDir, @(x) isnumeric(x) && isscalar(x) && (x~=0));
            parse(p,ftt,varargin{:});
            
            obj.dir = p.Results.dir;
            obj.ftt = round(ftt);
            
            d = length(obj.ftt.cores);
            obj.ys = cell(d, 1);
            if obj.dir > 0
                % last dimenion
                figeqk  = 1;
                for k = d:-1:1
                    nx  = obj.ftt.oneds{k}.num_nodes;
                    rkm = size(obj.ftt.cores{k}, 1);
                    rk  = size(obj.ftt.cores{k}, 3);
                    % push the previous dimenion into the current ftt by modifying the
                    % coefficients, then ys{k} is used to replace the ftt core at the
                    % k-th coordinate for obtaining the integrated function over the
                    % last d-k dimensions, ys{k} is a rk-1 by nx matrix
                    obj.ys{k} = reshape(reshape(obj.ftt.cores{k}, rkm*nx, rk)*figeqk, rkm, nx);
                    figeqk  = oned_integral(obj.ftt.oneds{k}, obj.ys{k}')';
                end
                obj.z = figeqk;
            else
                % first dimenion
                fileqk = 1;
                for k = 1:d
                    nx  = obj.ftt.oneds{k}.num_nodes;
                    rkm = size(obj.ftt.cores{k}, 1);
                    rk  = size(obj.ftt.cores{k}, 3);
                    % push the previous dimenion into the current ftt
                    % ys{k} is a nx by rk matrix
                    obj.ys{k} = reshape(fileqk*reshape(obj.ftt.cores{k}, rkm, nx*rk), nx, rk);
                    fileqk  = oned_integral(obj.ftt.oneds{k}, obj.ys{k});
                end
                obj.z = fileqk;
                % end
            end
            
            obj.oned_cdfs = cell(size(obj.ftt.oneds));
            for i = 1:length(ftt.oneds)
                obj.oned_cdfs{i} = CDFconstuctor(ftt.oneds{i}); % not implemented
            end
        end
        
        function [r,f] = eval_irt(obj, z)
            %Using IRT to sample from the target pdf approximated by firt
            %
            %%%%
            %Inputs:
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
            
            d = length(obj.ftt.cores);
            n = size(z,2);
            if size(z,1) ~= d
                error('Input variable must be d x n')
            end
            r = zeros(d,n);
            
            if obj.dir > 0
                frl = ones(n,1);
                for k = 1:d
                    rkm = size(obj.cores{k}, 1);
                    pk  = obj.ys{k}'*frl';
                    r(k,:)  = sample_oned_cdf(obj.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, r(k,:));
                    jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                    ii  = repmat((1:n)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                f   = frl'/obj.z;
            else
                frg = ones(1,n);
                for k = d:-1:1
                    rk  = size(obj.cores{k}, 3);
                    pk  = obj.ys{k}*frg;
                    r(k,:)  = sample_oned_cdf(obj.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, r(k,:));
                    ii  = reshape(1:rk*n, [], 1);
                    jj  = reshape(repmat(1:n, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*n, n);
                    frg = T2'*B;
                end
                f   = frg/obj.z;
            end
        end
        
        function [r,f] = eval_cond_irt(obj, x, z)
            %Using IRT to draw conditional samples from the target pdf approximated by
            %firt
            %
            %%%%
            %Inputs:
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
            
            
            d  = length(obj.cores);
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
            if obj.dir > 0
                %from the last dimension, marginalise to the first
                %order of the samples (x, r)
                %first evaluate the marginal for x
                
                frl = ones(nr,1);
                for k = 1:(dx-1)
                    rkm = size(obj.cores{k}, 1);
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, x(k,:));
                    jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
                    ii  = repmat((1:nx)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                
                k   = dx;
                rkm = size(obj.cores{k}, 1);
                % evaluate the updated basis function
                % for the marginal
                fxk = eval_oned_nodes2x(obj.oneds{k}, x(k,:), obj.ys{k}');
                % mutiply together
                fm  = sum(frl.*fxk, 2)';
                
                %for the frl, evaluate the updated basis function
                T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, x(k,:));
                jj  = reshape(reshape(1:rkm*nx, rkm, nx)', [], 1);
                ii  = repmat((1:nx)', 1, rkm);
                B   = sparse(ii(:), jj(:), frl(:), nx, rkm*nx);
                frl = B*T2;
                
                %then generate the conditional samples
                for j = 1:dr
                    k = dx+j;
                    rkm = size(obj.cores{k}, 1);
                    pk  = obj.ys{k}'*frl';
                    
                    r(j,:)  = sample_oned_cdf(obj.oned_cdfs{k}, pk, z(j,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, r(j,:));
                    jj  = reshape(reshape(1:rkm*nr, rkm, nr)', [], 1);
                    ii  = repmat((1:nr)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), nr, rkm*nr);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
                f   = frl';
            else
                %from the first dimension, marginalise to the last
                %order of the samples (r, x)
                %first evaluate the marginal for x
                
                frg = ones(1,nr);
                for j = dx:-1:2
                    k   = dr + j;
                    rk  = size(obj.cores{k}, 3);
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, x(j,:));
                    ii  = reshape(1:rk*nx, [], 1);
                    jj  = reshape(repmat(1:nx, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*nx, nx);
                    frg = T2'*B;
                end
                
                k   = dr+1;
                rk  = size(obj.cores{k}, 3);
                % evaluate the updated basis function
                
                % for marginal
                fxk = eval_oned_nodes2x(obj.oneds{k}, x(1,:), obj.ys{k});
                % mutiply together
                fm  = sum(fxk'.*frg, 1);
                
                % for frg
                T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, x(1,:));
                ii  = reshape(1:rk*nx, [], 1);
                jj  = reshape(repmat(1:nx, rk, 1), [], 1);
                B   = sparse(ii, jj, frg(:), rk*nx, nx);
                frg = T2'*B;
                
                %then generate the conditional samples
                for k = dr:-1:1
                    rk  = size(obj.cores{k}, 3);
                    pk  = obj.ys{k}*frg;
                    r(k,:)  = sample_oned_cdf(obj.oned_cdfs{k}, pk, z(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, r(k,:));
                    ii  = reshape(1:rk*nr, [], 1);
                    jj  = reshape(repmat(1:nr, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*nr, nr);
                    frg = T2'*B;
                end
                f   = frg;
            end
            f = f./fm;
            
        end
        
        function z = eval_rt(obj, r)
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
            
            
            d = length(obj.cores);
            n = size(r,2);
            if size(r,1) ~= d
                error('Input variable must be d x n')
            end
            z = zeros(d,n);
            
            if obj.dir > 0
                frl = ones(n,1);
                for k = 1:d
                    rkm = size(obj.cores{k}, 1);
                    pk  = obj.ys{k}'*frl';
                    z(k,:)  = eval_oned_cdf(obj.oned_cdfs{k}, pk, r(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_213(obj.oneds{k}, obj.cores{k}, r(k,:));
                    jj  = reshape(reshape(1:rkm*n, rkm, n)', [], 1);
                    ii  = repmat((1:n)', 1, rkm);
                    B   = sparse(ii(:), jj(:), frl(:), n, rkm*n);
                    frl = B*T2;
                    %frl(isnan(frl)) = 0;
                end
            else
                frg = ones(1,n);
                for k = d:-1:1
                    rk  = size(obj.cores{k}, 3);
                    pk  = obj.ys{k}*frg;
                    z(k,:)  = eval_oned_cdf(obj.oned_cdfs{k}, pk, r(k,:));
                    
                    % evaluate the updated basis function
                    T2  = eval_oned_core_231(obj.oneds{k}, obj.cores{k}, r(k,:));
                    ii  = reshape(1:rk*n, [], 1);
                    jj  = reshape(repmat(1:n, rk, 1), [], 1);
                    B   = sparse(ii, jj, frg(:), rk*n, n);
                    frg = T2'*B;
                end
            end
        end
        
    end
    
end