function [xcoord, ycoord, BlobsPerLayer] = Discretize_Single_Shell_Disk(a,s)
%Returns BlobsPerLayer and the coordinates of each blob for a single disk.
%   Detailed explanation goes here

NR = 1;

BlobsPerLayer = zeros([1, NR]);

BlobsPerLayer(1) = 200; 


Nblobs = sum(BlobsPerLayer);
%%% we need to find the coordinates of the blobs
xcoord = zeros([1, Nblobs]);
ycoord = zeros([1, Nblobs]);


for j=1:BlobsPerLayer(1)
    delta_phi = 2*pi/BlobsPerLayer(1);
    xcoord(j) = a * cos((j-1) * delta_phi); 
    ycoord(j) = a * sin((j-1) * delta_phi);
end
end
