function r = root_finding(func, tol, method, a, b)

switch method
    case{'Regula_Falsi'}
        r = regula_falsi(func, tol, a, b);
    case{'Illinois'}
        r = illinois(func, tol, a, b);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = regula_falsi(func, tol, a, b)


fa = func(a) ;
fb = func(b) ;

if sum(sign(fb.*fa) ~= -1)
    disp('Root finding: initial guesses on one side')
end

c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cold    = inf;

while ( norm(c-cold, inf) >= tol )
    cold    = c;
    fc  = func(c) ;
    
    I1  = (fc < 0);
    I2  = (fc > 0);
    I3  = ~I1 & ~I2;
    a   = I1.*c + I2.*a + I3.*c;
    b   = I1.*b + I2.*c + I3.*c;
    fa  = I1.*fc + I2.*fa + I3.*fc;
    fb  = I1.*fb + I2.*fc + I3.*fc;
    step    = -fb.*(b - a)./(fb - fa);
    step(isnan(step)) = 0;
    c = b + step;
    
    %norm(fc, inf)
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = illinois(func, tol, a, b)

fa = func(a) ;
fb = func(b) ;

c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cold    = inf;

side = zeros(size(a));

while ( norm(c - cold, inf) >= tol )
    cold    = c;
    fc  = func(c) ;
    
    I1  = (fc < 0);
    I2  = (fc > 0);
    I3  = ~I1 & ~I2;
    a   = I1.*c + I2.*a + I3.*c;
    b   = I1.*b + I2.*c + I3.*c;
    fa  = I1.*fc + I2.*fa + I3.*fc;
    fb  = I1.*fb + I2.*fc + I3.*fc;
    
    fb(side(I1) == -1) = fb(side(I1) == -1)/2;
    side(I1) = -1;
    fa(side(I1) == 1) = fa(side(I1) == 1)/2;
    side(I2) = 1;
    
    step = -fb.*(b - a)./(fb - fa);
    step(isnan(step)) = 0;
    c = b + step;
end

end

