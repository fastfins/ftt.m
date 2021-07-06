function c = regula_falsi(func, tol, a, b)


fa = func(a) ;
fb = func(b) ;

if sum(sign(fb.*fa) ~= -1)
    disp('Root finding: initial guesses on one side')
end

c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
cold    = inf;

i = 2;
while ( norm(c-cold, Inf) >= tol )
    cold    = c;
    fc  = func(c);
    if norm(fc, Inf) < tol
        break;
    end
    
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
    i = i+1;
end
disp(i)

end