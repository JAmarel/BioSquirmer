tic

%Simulation
T = 20;
dt = 1;


%Discretization
a = 1;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.1 * a;          %%% radial spacing between neighboring blobs
epsilon = s/8;       %%% radius of blobs

%Enclosure
R = 5*a;   %%%Radius of enclosure
d = 2*s;    %%%Circumferential Enclosure Blob Spacing

%Initial Conditions
r_o = 4*a;          %%% Radial coordinate of beast cm from center of enclosure
phi_o = 0*pi/2;     %%%Angle coordinate of beast cm from center of enclosure
theta_o = 2*pi/4;   %%% Beast intial orientation (head direction)

%Coordinates of beast blobs in beast frame.
%Beast frame has its head on the x axis at y = 0.
[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

Nblobs = sum(BlobsPerLayer); %%% Number of blobs in the beast
NRim = BlobsPerLayer(end);   %%% Number of blobs in the outermost beast layer

%Enclosure Blob Coordinates from the origin in the enclosure frame.
[x_Enc, y_Enc] = DiscretizeEnclosure(R,d); 


%Prescribe wave in the enclosure frame, taking into account theta_o rotation of the head.
[VxRim, VyRim, B1] = PrescribeWave_Orient(NRim,theta_o);

%Translate beast coordinates with (r_o, phi_o) and head facing direction (theta_o)
%Into the enclosure frame
[xcoord, ycoord, x_head, y_head] = Orient_disk(xcoord, ycoord, r_o, phi_o, theta_o, BlobsPerLayer);
%xcoord and ycoord now label beast blob coordinates in the enclosure frame.


[Ux_history, Uy_history, W_history, x_history, y_history, theta_history, Angles, x_cm_history, y_cm_history] = ...
    TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, r_o, phi_o);


toc

%%Plotting
figure(1)
plot(x_cm_history(1), y_cm_history(1), 'ko','LineWidth', 1) %Begin at black
hold on
plot(x_cm_history(end), y_cm_history(end), 'ro','LineWidth', 1) %End at red
plot(x_cm_history(2:end-1), y_cm_history(2:end-1), '.')
daspect([1,1,1])
%plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'r.','LineWidth', 1)
plot(x_Enc, y_Enc, 'go', 'LineWidth', 1)
%plot(x_head, y_head, 'ko', 'LineWidth', 1)
axis off
hold off