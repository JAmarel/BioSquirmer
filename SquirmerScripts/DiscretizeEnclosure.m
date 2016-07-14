function [xcoord, ycoord] = DiscretizeEnclosure(R,d)
%Returns the coordinates of each blob for a single layer enclosure.
%   R is the radius of the enclosure
%   d is the circumferential blob spacing

Nblobs = round(2*pi*R/d,0); %Packs as many blobs in as possible according to R and d.
%%% we need to find the coordinates of the blobs
xcoord = zeros([1, Nblobs]);
ycoord = xcoord;

xcoord(1)=0; %%% blob in the center of the circle
ycoord(1)=0; %%% blob in the center of the circle

delta_phi = 2*pi/Nblobs;
        
for i=1:Nblobs
    xcoord(i) =  R * cos((i-1) * delta_phi); 
    ycoord(i) =  R * sin((i-1) * delta_phi);
end
end

