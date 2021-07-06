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
        fx = sum((fxl*obj.ms{dz+1}).^2, 2)'/obj.z;
    else
        fx = fxl'.^2/obj.z;
    end
else
    % marginalised from the left
    fxg = eval_block(obj, x, obj.int_dir);
    if dz < d
        ie = (d-dz)+1;
        fx = sum((obj.ms{ie-1}*fxg).^2, 1)/obj.z;
    else
        fx = fxg.^2/obj.z;
    end
end

end
