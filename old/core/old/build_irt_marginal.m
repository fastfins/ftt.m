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