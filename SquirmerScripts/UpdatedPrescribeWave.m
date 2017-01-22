function [VxRim, VyRim] = UpdatedPrescribeWave(NRim, B1, B2, theta_o)
%Prescribes single squirmer wave when the beast head is facing theta_o
%where theta_o is the angle of the beasts face as measured from the
%horizontal x axis of the ring.

%Also nondimensionalizes the wave by dividing through by B1/2 (or B2/2 if B1=0)


VxRim     = zeros([NRim, 1]);     %%% x-component of tangential velocity at the rim
VyRim     = zeros([NRim, 1]);     %%% y-component of tangential velocity at the rim

for i=1:NRim
        
        angle = (i-1) * 2 * pi/NRim;
        VRimTheta = B1 * sin(angle - theta_o) + B2 * sin(2 * (angle - theta_o));
        VxRim(i) = -VRimTheta * sin(angle); 
        VyRim(i) = VRimTheta * cos(angle);
end

%Nondimensionalize
if B1 == 0
    VxRim = VxRim./abs((B2/2));
    VyRim = VyRim./abs((B2/2));
elseif B1==0 && B2==0
    VxRim = VxRim;
    VyRim = VyRim;
else
    VxRim = VxRim./abs((B1/2));
    VyRim = VyRim./abs((B1/2));
end

end

