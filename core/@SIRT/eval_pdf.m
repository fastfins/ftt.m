function fx = eval_pdf(obj, x)
% Evaluate the marginalise pdf represented by ftt
%
%   x   - input variables
%
%   fx  - marginal density at x

dz = size(x,1);
d = length(obj.cores);
if obj.marginal_direction > 0
    % marginalised from the right
    fxl = eval_block(obj, x, obj.marginal_direction);
    if dz < d
        fx = sum((fxl*obj.ms{dz+1}).^2, 2)'/obj.z;
    else
        fx = fxl'.^2/obj.z;
    end
else
    % marginalised from the left
    fxg = eval_block(obj, x, obj.marginal_direction);
    if dz < d
        ie = (d-dz)+1;
        fx = sum((obj.ms{ie-1}*fxg).^2, 1)/obj.z;
    else
        fx = fxg.^2/obj.z;
    end
end

end
