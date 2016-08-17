% Second inactive attempt. Pick a U and solve for F by summing the little f's.
% Active and Inactive disks must have velocity in the same direction (x)
tic
n = 20;
Radii = logspace(-1.5,3.5,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

ActiveSpeeds = zeros([1, n]);
RecipSpeeds = zeros([1, n]);

for i=1:n
    a = Radii(i);
    s = 0.1 * a;          %%% spacing between neighboring blobs
    epsilon = s/8;       %%% radius of the blob

    %Velocity of inactive disk
    Ux1 = 1; 
    Uy1 = 0;

    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

    Nblobs = sum(BlobsPerLayer); %%% total number of blobs 

    NR = length(BlobsPerLayer); %%% Number of radial layers
    NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

    [VxRim, VyRim, B1] = PrescribeWave(NRim);

    %Nondimensionalize
    VxRim = VxRim/(B1/2);
    VyRim = VyRim/(B1/2);

    %% Inactive
    [fx1, fy1, Matrix1] = solve_U_disk_inactive(xcoord, ycoord, epsilon, NRim, Ux1, Uy1);

    FxRim1 = fx1(end-NRim+1:end);
    FyRim1 = fy1(end-NRim+1:end);

    FxNet1 = sum(fx1); %%% x-component of net force on squirmer
    FyNet1 = sum(fy1); %%% y-component of net force on squirmer

    TorqueNet1 = dot(xcoord,fy1.') - dot(ycoord,fx1.');

    %% Active
    [fx2, fy2, Ux2, Uy2, W2, Matrix2] = solve_U_disk_rot(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);

    FxRim2 = fx2(end-NRim+1:end);
    FyRim2 = fy2(end-NRim+1:end);

    FxNet2 = sum(fx2); %%% x-component of net force on squirmer
    FyNet2 = sum(fy2); %%% y-component of net force on squirmer

    TorqueNet2 = dot(xcoord,fy2.') - dot(ycoord,fx2.');
    
    ActiveSpeeds(i) = Ux2;
    %% Recip Thm

    % hat corresponds to the conjugate scenario
    F_hat_x = Ux1/HPW_mobility(a);
    F_hat_y = Uy1/HPW_mobility(a);


    % These terms are from equation 12 in our notes
    LHS = F_hat_x*Ux2 + F_hat_y*Uy2;
    RHS = -(sum(VxRim.*FxRim1) + sum(VyRim.*FyRim1));
    Ratio = LHS/RHS; %Should be close to 1

    %This is true iff velocity is only in x
    Ux = (1/(FxNet1))*RHS; %Ux of the active swimmer as predicted by the recip thm
    RecipSpeeds(i) = Ux;

end
toc


%% Plotting
figure(1)
semilogx(Radii, ActiveSpeeds, 'k-')
hold on
semilogx(Radii, RecipSpeeds, 'bo')
hold off
title('Swimming Speed vs. Radius','FontSize',16,'FontWeight','bold')
xlabel('Log Scale Nondimensional Radius (a/l_s)')
ylabel('Nondimensional Swimming Speed U/(B1/2)')
