function [xcoord, ycoord, x_head, y_head] = Orient_disk(xcoord, ycoord, r_o,phi_o, theta_o, BlobsPerLayer)
%%Place beast somewhere (r,phi) with head facing direction (theta_o).

x_o = r_o*cos(phi_o); %%%Beast CM initial x position as seen in enclosure frame.
y_o = r_o*sin(phi_o); %%%Beast CM Initial y position


%Find the angles of beast blobs, as seen in the beast frame.
Angles = zeros([1, length(xcoord)]);
%Get Angles from arctan. [-pi/2,pi/2]. Break up quadrants.
for i=2:length(xcoord)
    if xcoord(i)>=0 && ycoord(i)>=0
        Angles(i) = atan(abs(ycoord(i)/xcoord(i)));
    elseif xcoord(i)<0 && ycoord(i)>=0
        Angles(i) = pi - atan(abs(ycoord(i)/xcoord(i)));
    elseif xcoord(i)<0 && ycoord(i)<0
        Angles(i) = pi + atan(abs(ycoord(i)/xcoord(i)));
    elseif xcoord(i)>=0 && ycoord(i)<0
        Angles(i) = 2*pi - atan(abs(ycoord(i)/xcoord(i)));
    end
end


Angles = Angles + theta_o;
Angles(1) = 0;

%Blob distances from beast center
BlobDistance = sqrt(xcoord.^2 + ycoord.^2); %Should be shell radii

%Rotate all coordinates according to the new angle
xcoord = BlobDistance.*cos(Angles);
ycoord = BlobDistance.*sin(Angles);
xcoord(1) = 0;
ycoord(1) = 0;

%Translate to enclosure frame
xcoord = xcoord + x_o;
ycoord = ycoord + y_o;

%Find head coordinates. I don't like the way this needs bpl.
%Nblobs = length(BlobsPerLayer); %%% total number of beast blobs 
NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

x_head = xcoord(end - NRim + 1);
y_head = ycoord(end - NRim + 1);

end

