function [xcoord, ycoord] = DiscretizeWall(R,s,h)
%Returns the coordinates of each blob for a single vertical layer wall of height h that is symmetrically placed a distance R from the origin.
%   R is the radius of the enclosure
%   s is the blob spacing

Nblobs = round(h/s,0); %Packs as many blobs in as possible according to h and s.


xcoord = R*ones([1, Nblobs]);
ycoord = linspace(-h/2,h/2,Nblobs);
end

