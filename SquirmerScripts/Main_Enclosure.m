tic

a = 10;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.1 * a;          %%% spacing between neighboring blobs
epsilon = s/8;       %%% radius of the blob

R = 10*a; %%%Radius of enclosure
d = 5*s;    %%%Circumferential Enclosure Blob Spacing

r = 3*a; %%% Radial coordinate of beast from center of enclosure
phi = pi/2; %%%Angle coordinate of beast from center of enclosure

theta_o = pi/4;   %%% Beast intial orientation
x_o = r*cos(phi); %%%Beast initial x position
y_o = r*sin(phi); %%%Initial y position

%Blob coordinates from beast center
[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

%Get Angles from arctan
BlobDistance = sqrt(xcoord.^2 + ycoord.^2); %Blob distances from beast center


xcoord = xcoord + x_o;
ycoord = ycoord + y_o;

[x_Enc, y_Enc] = DiscretizeEnclosure(R,d); %Enclosure Blob Coordinates

Nblobs = sum(BlobsPerLayer); %%% total number of beast blobs 

NR = length(BlobsPerLayer); %%% Number of radial layers
NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

x_head = xcoord(end - NRim + 1); %Check these to be sure
y_head = ycoord(end - NRim + 1);

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