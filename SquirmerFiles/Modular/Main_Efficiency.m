n = 20;
Radii = logspace(-2,2,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

Efficiencies = zeros([1, n]);

for i=1:n
    a = Radii(i);
    s= 0.1 * a;          %%% spacing between neighboring blobs
    epsilon = s/8;       %%% radius of the blob
    
    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);
    Nblobs = sum(BlobsPerLayer); %%% total number of blobs 
    
    NR = length(BlobsPerLayer); %%% Number of radial layers
    NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer
    
    [VxRim, VyRim, B1] = PrescribeWave(NRim);
    [fx, fy, Ux, Uy] = solve_U_disk(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);

    FxRim = fx(end-NRim+1:end);
    FyRim = fy(end-NRim+1:end);

    FxNet = sum(fx); %%% x-component of net force on squirmer
    FyNet = sum(fy); %%% y-component of net force on squirmer

    speed = Ux/(B1/2);   %%%% swimming velocity non-dimensionalized by B1/2

    efficiency = CalcEfficiency(FxRim, FyRim, VxRim, VyRim, a, speed);
    Efficiencies(i) = efficiency;
end

%%Plot the efficiency vs nondimensional radius
figure(1)
plot(Radii, Efficiencies, 'o')

figure(2)
semilogx(Radii, Efficiencies, 'o')
