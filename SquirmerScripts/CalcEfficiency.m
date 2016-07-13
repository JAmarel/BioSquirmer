function efficiency = CalcEfficiency(FxRim, FyRim, VxRim, VyRim, a, speed)
%Ppull/Pswim
%   Detailed explanation goes here
Pswim = dot(FxRim,VxRim) + dot(FyRim,VyRim);

mu = HPW_mobility(a);

Ppull = (speed^2)/mu;

efficiency = Ppull/Pswim;

end

