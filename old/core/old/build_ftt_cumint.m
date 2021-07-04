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
