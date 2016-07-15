tic

a = 10;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.1 * a;          %%% radial spacing between neighboring blobs
epsilon = s/8;       %%% radius of blobs

R = 10*a;   %%%Radius of enclosure
d = 5*s;    %%%Circumferential Enclosure Blob Spacing

r_o = 3*a;          %%% Radial coordinate of beast cm from center of enclosure
phi_o = 3*pi/2;     %%%Angle coordinate of beast cm from center of enclosure
theta_o = 3*pi/4;   %%% Beast intial orientation (head direction)

%Blob coordinates from beast center
[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

Nblobs = sum(BlobsPerLayer);
NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

[VxRim, VyRim, B1] = PrescribeWave(NRim);

%Place beast somewhere in enclosure (r,phi) with head facing direction (theta_o).
[xcoord, ycoord, x_head, y_head] = Orient_disk(xcoord, ycoord, r_o,phi_o, theta_o, BlobsPerLayer);

[x_Enc, y_Enc] = DiscretizeEnclosure(R,d); %Enclosure Blob Coordinates

[fx, fy, Ux, Uy, W, Ux_Enc, Uy_Enc, Matrix] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim);


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