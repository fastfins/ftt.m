function fx = eval_pdf(obj, x)
% Evaluate the normalised (marginal) pdf represented by squared FTT.
%   f = EVAL_PDF(irt, x)
%
%   x - input variables
%   f - marginal density at x

dz = size(x,1);
d = length(obj.cores);
if obj.int_dir > 0
    % marginalised from the right
    fxl = eval_block(obj, x, obj.int_dir);
    if dz < d
        fx = sum((fxl*obj.ms{dz+1}).^2, 2)';
    else
        fx = fxl'.^2;
    end
    fx_ref = ones(1,size(x,2));
    for k = 1:dz 
        fk_ref = eval_pdf(obj.oned_refs{k}, x(k,:));
        fx_ref = fx_ref.*fk_ref;
    end
else
    % marginalised from the left
    fxg = eval_block(obj, x, obj.int_dir);
    if dz < d
        ie = (d-dz)+1;
        fx = sum((obj.ms{ie-1}*fxg).^2, 1);
    else
        fx = fxg.^2;
    end
    fx_ref = ones(1,size(x,2));
    for j = 1:dz
        k = j + (d-dz);
        fj_ref = eval_pdf(obj.oned_refs{k}, x(j,:));
        fx_ref = fx_ref.*fj_ref;
    end
end
fx = (fx+obj.tau*fx_ref)/obj.z;

end
