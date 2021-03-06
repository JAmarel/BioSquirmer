function [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s)
%Returns BlobsPerLayer and the coordinates of each blob for a single disk.
NR = floor(a/s) + 1;  %%% number of layers along radial direction (+1 b/c of one blob that sits in the center)

%%% floor rounds the number to the nearest integer less or equal to that number
%%% If a/s is near a half integer value, floor has trouble rounding
%%% consitently, use round instead. And vcv if the value is near whole
%%% integer.



BlobsPerLayer = zeros([1, NR]);

BlobsPerLayer(1) = 1; %%% one blob in the center of the circle

for i = 2: NR
    BlobsPerLayer(i) = round( (2 * pi * (i - 1))); 
end


Nblobs = sum(BlobsPerLayer);
%%% we need to find the coordinates of the blobs
xcoord = zeros([1, Nblobs]);
ycoord = zeros([1, Nblobs]);

xcoord(1)=0; %%% blob in the center of the circle
ycoord(1)=0; %%% blob in the center of the circle

index = 2; % index equivalent to i?
for i=2:NR
    for j=1:BlobsPerLayer(i)
        delta_phi = 2*pi/BlobsPerLayer(i);
        xcoord(index) = (i-1) * s * cos((j-1) * delta_phi); 
        ycoord(index) = (i-1) * s * sin((j-1) * delta_phi);
        index = index + 1;
    end
end
end

