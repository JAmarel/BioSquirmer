function [VxRim, VyRim, B1] = PrescribeWave_Orient(NRim,theta_o)
%Prescribes single squirmer wave in the enclosure frame.
%   This allows for some intial beast head orientation

B1 = 0.1;            %%% tangential velocity strength (streaming)
B2 = 0;   

VRimTheta = zeros([NRim, 1]);     %%% tangential velocity at the rim of the disk
VxRim     = zeros([NRim, 1]);     %%% x-component of tangential velocity at the rim
VyRim     = zeros([NRim, 1]);     %%% y-component of tangential velocity at the rim
Angles    = zeros([NRim,1]);      %%% Angles corresponding to rim blobs.

%First entry corresponds to the head blob.

for i=1:NRim
    
    angle = (i-1) * 2 * pi/NRim;
    VRimTheta(i) = B1*sin(angle) + B2*sin(2 * angle);
    
    VxRim(i)  = -VRimTheta(i) * sin(angle); 
    VyRim(i)  = VRimTheta(i) * cos(angle);
    
    Angles(i) = angle;
end

%Now Rotate velocities according to theta_o into lab frame.
VxRimNew = VxRim*cos(theta_o) - VyRim*sin(theta_o);
VyRimNew = VxRim*sin(theta_o) + VyRim*cos(theta_o);

VxRim = VxRimNew;
VyRim = VyRimNew;



end

