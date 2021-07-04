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

