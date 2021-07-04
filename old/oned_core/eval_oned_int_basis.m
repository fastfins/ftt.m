function b = eval_oned_int_basis(type, domain, order, x)
%Evaluate the integral of chebyshev basis functions at input x
%
%Tiangang Cui, August, 2019

if order < 1
    disp('order is less than two, this is not needed')
    return
end

[x,x2z] = domain2reference(x(:), domain);

switch type
    case{'Chebyshev1st'}
        theta = real(acos(x));
        normalising = reshape([1,sqrt(2)*ones(1,order)]/sqrt(pi), 1,[]);
        %
        b = zeros(length(x),order+1);
        b(:,1) = x; % n = 0
        b(:,2) = 0.5*x.^2; % n = 1
        if order > 2
            % int(T_n) = T_(n+1) ( n/(n^2-1) ) - x T_n / (n-1), n > 1
            n = 2:order;
            % evaluate two temporary cheby polys
            tmp1 = cos( theta.*(n+1) ).*(n./(n.^2-1));
            tmp2 = cos( theta.*n ).*(x./(n-1));
            b(:,n+1) = (tmp1-tmp2);
        end
        b = b.*normalising;
    case{'Chebyshev2nd'}
        theta = real(acos(x));
        % the normalising constant used
        normalising = sqrt(2/pi);
        % int(U_n) = T_(n+1) / (n+1), n = 0, ..., order
        n = 0:order;
        b = cos( theta.*(n+1) ) .* (normalising./(n+1));
        
    case{'Fourier'}
        
        m   = order+1;
        n   = m*2;
        is  = 2:m;
        ic  = (m+1):(n-1);
        
        b = zeros(length(x), n);
        b(:,1)  = x(:)*sqrt(0.5);
        b(:,is) = - cos( x(:)*(1:order)*pi )./ (pi*(1:order));
        b(:,ic) = sin( x(:)*(1:order)*pi )  ./ (pi*(1:order));
        b(:,n)  = sin( x(:)*(m*pi) )*sqrt(0.5)/(pi*m);
        
end

b = b.*x2z(:);

end
% use the following normalisation
%
%
%
%normalising2nd = sqrt(2/pi);