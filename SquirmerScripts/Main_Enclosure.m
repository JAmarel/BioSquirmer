tic

%Simulation
T = 2;
dt = .5;


%Discretization
a = 1;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.1 * a;          %%% radial spacing between neighboring blobs
epsilon = s/8;       %%% radius of blobs

%Enclosure
R = 5*a;   %%%Radius of enclosure
d = 5*s;    %%%Circumferential Enclosure Blob Spacing

%Initial Conditions
r_o = 1*a;          %%% Radial coordinate of beast cm from center of enclosure
phi_o = 2*pi/2;     %%%Angle coordinate of beast cm from center of enclosure
theta_o = 3*pi/4;   %%% Beast intial orientation (head direction)

%Blob coordinates from beast center
[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

%Enclosure Blob Coordinates. Enclosure is center at the origin.
[x_Enc, y_Enc] = DiscretizeEnclosure(R,d); 


Nblobs = sum(BlobsPerLayer); %%% number of blobs in the beast
NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost beast layer


[VxRim, VyRim, B1] = PrescribeWave(NRim); %This needs to be recalculated to account for theta_o and dt's.

%Place beast somewhere in enclosure (r_o,phi_o) with head facing direction (theta_o).
[xcoord, ycoord, x_head, y_head] = Orient_disk(xcoord, ycoord, r_o, phi_o, theta_o, BlobsPerLayer);


[Ux_history, Uy_history, W_history, x_history, y_history, theta_history, Angles] = ...
    TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, r_o, phi_o);


toc

%Plot the position of all blobs
figure(1)
plot(xcoord, ycoord, '.')
daspect([1,1,1])
hold on
plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'r.', 'LineWidth', 1)
plot(x_Enc, y_Enc, 'go', 'LineWidth', 1)
plot(x_head, y_head, 'ko', 'LineWidth', 1)
axis off
hold off