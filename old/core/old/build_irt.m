function ftt_irt = build_irt(ftt, varargin)
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

ftt_irt.dir     = dir;
ftt_irt.cores   = ftt.cores;
ftt_irt.oneds   = ftt.oneds;
ftt_irt.ng_flag = ftt.ng_flag;
[ftt_irt.z, ftt_irt.ys]  = build_ftt_cumint(ftt, dir);

ftt_irt.oned_cdfs = cell(size(ftt.oneds));
for i = 1:length(ftt.oneds)
    ftt_irt.oned_cdfs{i} = setup_oned_cdf(ftt.oneds{i}, 'squared', ftt.ng_flag);
end


end