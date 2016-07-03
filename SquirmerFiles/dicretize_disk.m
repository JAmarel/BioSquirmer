function dicretize_disk

a = 1;  %%% radius of the disk
s= 0.1 *a; %%% spacing between neighboring blobs
NR = floor(a/s) + 1; %%% number of layers along radial direction

%%% floor rounds the number to the nearest integer less or equal to that
%%% number; one blob will be in the center of the circle

BlobsPerLayer = zeros([1, length(NR)]);

BlobsPerLayer(1) = 1; %%% one blob in the center of the circle
Nblobs = 1;
for i = 2: NR
    BlobsPerLayer(i) = round( 2 * pi * (i - 1) ); 
    Nblobs = Nblobs + BlobsPerLayer(i);
end

Nblobs %%% total number of blobs 
BlobsPerLayer(NR)  %%% number of blobs in the outermost layer

%%% we need to find the coordinates of the blobs
xcoord = zeros([1, Nblobs]);
ycoord = xcoord;

xcoord(1)=0; %%% blob in the center of the circle
ycoord(1)=0; %%% blob in the center of the circle


index = 2;
for i=2:NR
    for j=1:BlobsPerLayer(i)
        delta_phi = 2*pi/BlobsPerLayer(i);
        xcoord(index) = (i-1) * s * cos( (j-1) * delta_phi); 
        ycoord(index) = (i-1) * s * sin( (j-1) * delta_phi);
        index = index + 1;
    end
end

figure(1)
plot(xcoord, ycoord, 'o')
daspect([1,1,1])
hold on
plot(xcoord(Nblobs-BlobsPerLayer(NR)+1:end), ycoord(Nblobs-BlobsPerLayer(NR)+1:end), 'ro')
hold off




end