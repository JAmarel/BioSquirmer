%%Calculates/plots the speed for various nondimensional beast radii. Unbounded
n = 20; %Number of radii to sample
Radii = logspace(-1.5,3.5,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.
B1 = 1;
B2 = 0;

Speeds = zeros([1, n]);

for i=1:n
    a = Radii(i);
    s = 0.08 * a;          %%% spacing between neighboring blobs
    epsilon = s/8;       %%% radius of the blob
    
    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);
    
    Nblobs = sum(BlobsPerLayer); %%% total number of blobs 
    NR = length(BlobsPerLayer); %%% Number of radial layers
    NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer
    
    [VxRim, VyRim] = PrescribeWave(NRim, B1, B2);
    
    %Nondimensionalize
    VxRim = VxRim/(B1/2);
    VyRim = VyRim/(B1/2);
    
    
    [fx, fy, Ux, Uy] = solve_U_disk(xcoord, ycoord, epsilon, VxRim, VyRim, NRim); %Currently not calling the rotating solver.
    %Use solve_U_disk_rot to allow rotational freedom. They agree when
    %unbounded


    speed = sqrt(Ux^2 + Uy^2);

    Speeds(i) = speed;
end

figure(1)
semilogx(Radii, Speeds, 'o')
title('Swimming Speed vs. Radius','FontSize',16,'FontWeight','bold')
xlabel('Nondimensional Beast Radius [a/l_s]')
ylabel('Nondimensional Swimming Speed [v/(B_1/2)] ')
saveas(gcf,'SpeedvsRadius.png')
saveas(gcf,'SpeedvsRadius.eps')

