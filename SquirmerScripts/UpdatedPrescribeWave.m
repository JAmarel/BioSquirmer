function [VxRim, VyRim] = UpdatedPrescribeWave(NRim, B1, B2, theta_o)
%Prescribes single squirmer wave when the beast head is facing theta_o
%according to the enclosure frame.


VxRim     = zeros([NRim, 1]);     %%% x-component of tangential velocity at the rim
VyRim     = zeros([NRim, 1]);     %%% y-component of tangential velocity at the rim

for i=1:NRim
        
        angle = (i-1) * 2 * pi/NRim;
        VRimTheta = 2 * sin(angle - theta_o) + 2 * (B2/B1) * sin(2 * (angle - theta_o)); %%% Non-dimensionalized by B1/2
        VxRim(i) = -VRimTheta * sin(angle); 
        VyRim(i) = VRimTheta * cos(angle);
end

end

