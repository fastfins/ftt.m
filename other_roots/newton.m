function c = newton(func, tol, a)

cold = a;
[f,df] = func(cold);
c = cold - f./df;

i = 1;
while ( norm(c-cold, Inf) >= tol )
    cold = c;
    [f,df] = func(cold);
    if norm(f, Inf)<tol
        break;
    end
    c = cold - f./df;
    %norm(f, inf)
    i = i+1;
end
disp(i)

end