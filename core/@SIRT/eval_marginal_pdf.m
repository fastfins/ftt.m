function fx = eval_marginal_pdf(obj, x)
% Evaluate the marginalise pdf represented by ftt
%
%   x   - input variables
%
%   fx  - marginal density at x

if obj.marginal_direction > 0
    % marginalised from the right
    fxl = eval_block(obj, x, obj.marginal_direction);
    if size(x,1) < length(obj.cores)
        fx = sum((fxl*obj.ms{dz}).^2, 2)'/obj.z;
    else
        fx = fxl'.^2/obj.z;
    end
else
    % marginalised from the left
    fxg = eval_block(obj, x, obj.marginal_direction);
    if size(x,1) < length(obj.cores)
        ie = (d-dz)+1;
        fx = sum((obj.ms{ie}*fxg).^2, 1)/obj.z;
    else
        fx = fxg.^2/obj.z;
    end
end

end
