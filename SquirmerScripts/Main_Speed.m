%%Calculates/plots the speed for various beast radii without rotation
n = 10;
Radii = logspace(-1.5,3.5,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

Speeds = zeros([1, n]);

for i=1:n
    a = Radii(i);
    s = 0.08 * a;          %%% spacing between neighboring blobs
    epsilon = s/8;       %%% radius of the blob
    
    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);
    
    Nblobs = sum(BlobsPerLayer); %%% total number of blobs 
    NR = length(BlobsPerLayer); %%% Number of radial layers
    NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer
    
    [VxRim, VyRim, B1] = PrescribeWave(NRim);
    
    
    [fx, fy, Ux, Uy] = solve_U_disk(xcoord, ycoord, epsilon, VxRim, VyRim, NRim); %Currently not calling the rotating solver.
    
    %Nondimensionalizing.
    Ux = Ux/(B1/2);
    Uy = Uy/(B1/2);


    speed = sqrt(Ux^2 + Uy^2);

    Speeds(i) = speed;
end

figure(2)
semilogx(Radii, Speeds, 'o')
title('Swimming Speed vs. Radius','FontSize',16,'FontWeight','bold')
xlabel('Log Scale Nondimensional Radius (a/l_s)')
ylabel('Nondimensional Swimming Speed')
%saveas(gcf,'SpeedvsRadius.png')

