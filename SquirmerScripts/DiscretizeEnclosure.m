function [xcoord, ycoord] = DiscretizeEnclosure(R,d)
%Returns the coordinates of each blob for a single layer enclosure.
%   R is the radius of the enclosure
%   d is the circumferential blob spacing
Nblobs = round(2*pi*R/d,0); %Packs as many blobs in as possible according to R and d.


xcoord = zeros([1, Nblobs]);
ycoord = zeros([1, Nblobs]);

delta_phi = 2*pi/Nblobs;

% Find the coordinates of the blobs  
for i=1:Nblobs
    xcoord(i) =  R * cos((i-1) * delta_phi); %(i-1) to begin at angle 0.
    ycoord(i) =  R * sin((i-1) * delta_phi);
end

% Adding another enclosure layer radially separated (outward) by 2d 
xcoord = [xcoord (1+d)*xcoord (1+2*d)*xcoord];
ycoord = [ycoord (1+d)*ycoord (1+2*d)*ycoord];

end

