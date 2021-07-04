function locx = process_bnd(dif, gs)
%Process the boudnary to avoid divide by zero
%
%Tiangang Cui, August, 2019

if gs < eps
    locx = 0;
else
    locx = dif/gs;
end

end