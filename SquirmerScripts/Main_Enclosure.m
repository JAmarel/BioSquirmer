tic

%Simulation
T = 400;
dt = .5;


%Discretization
a = 1;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.1 * a;          %%% radial spacing between neighboring blobs
epsilon = s/8;       %%% radius of blobs

%Enclosure
R = 10*a;   %%%Radius of enclosure
d = 2*s;    %%%Circumferential Enclosure Blob Spacing

%Initial Conditions
r_o = -8*a;          %%% Radial coordinate of beast cm from center of enclosure
phi_o = 0*pi/4;     %%%Angle coordinate of beast cm from center of enclosure
theta_o = pi/4;   %%% Beast intial orientation (head direction)

%Coordinates of beast blobs in beast frame.
%Beast frame has its head on its x axis at y = 0.
[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

%Rotate coordinates according to theta_o
xcoordNew = xcoord*cos(theta_o) - ycoord*sin(theta_o);
ycoordNew = xcoord*sin(theta_o) + ycoord*cos(theta_o);
xcoord = xcoordNew;
ycoord = ycoordNew;

Nblobs = sum(BlobsPerLayer); %%% Number of blobs in the beast
NRim = BlobsPerLayer(end);   %%% Number of blobs in the outermost beast layer

%Enclosure Blob Coordinates from the origin in the enclosure frame.
[x_Enc, y_Enc] = DiscretizeEnclosure(R,d); 


%Prescribe wave in the beast frame.
[VxRim, VyRim, B1] = PrescribeWave(NRim);

%Now Rotate velocities according to theta_o into lab frame.
VxRimNew = VxRim*cos(theta_o) - VyRim*sin(theta_o);
VyRimNew = VxRim*sin(theta_o) + VyRim*cos(theta_o);
VxRim = VxRimNew;
VyRim = VyRimNew;

%Translate beast CM to (r_o, phi_o) in enclosure frame
x_o = r_o*cos(phi_o); %%%Beast CM initial x position as seen in enclosure frame.
y_o = r_o*sin(phi_o); %%%Beast CM Initial y position

xcoord = xcoord + x_o;
ycoord = ycoord + y_o;

%Grab the head coordinate for plotting purposes
x_head = xcoord(end - NRim + 1);
y_head = ycoord(end - NRim + 1);



[Ux_history, Uy_history, W_history, x_history, y_history, theta_history, Angles, x_cm_history, y_cm_history] = ...
    TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, r_o, phi_o);


toc
%%Plotting
figure(1)
plot(x_cm_history(1), y_cm_history(1), 'Marker', 'o', 'MarkerSize', 2, ...
     'MarkerFaceColor', 'green'); %Begin at green
hold on
plot(x_cm_history(end), y_cm_history(end), 'Marker', 'o', 'MarkerSize', 2, ...
     'MarkerFaceColor', 'red'); %End at red
plot(x_cm_history(2:end-1), y_cm_history(2:end-1),'Marker', '.', 'MarkerSize', 2, ...
     'MarkerFaceColor', 'black');
daspect([1,1,1])
%plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'r.','LineWidth', 1)
plot(x_Enc, y_Enc, 'Marker', '.', 'MarkerSize', 1, ...
     'MarkerFaceColor', 'black');
%plot(x_head, y_head, 'ko', 'LineWidth', 1)
axis off
hold off