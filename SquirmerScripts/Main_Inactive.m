% Second inactive attempt. Pick a U and solve for F by summing the little f's.
% Active and Inactive disks must have velocity in the same direction (x)
tic
n = 20;
Radii = logspace(-1.5,3.5,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

ActiveSpeeds = zeros([1, n]);
RecipSpeeds = zeros([1, n]);

for i=1:n
    a = Radii(i);
    s = 0.1 * a;         %%% spacing between neighboring blobs
    epsilon = s/8;      %%% radius of the blob
    theta_o = 0;       %%% Beast intial orientation (head direction)

    %Velocity of inactive disk
    U1 = 1;
    Ux1 = U1*cos(theta_o); 
    Uy1 = U1*sin(theta_o);

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
    F_hat = U1/HPW_mobility(a);
    
    %FNet1 gives perfect agreement, but is incorrect?
    FNet1 = sqrt(FxNet1^2 + FyNet1^2);


    % These terms are from equation 12 in our notes
    U = (-1/(F_hat))*(sum(VxRim.*FxRim1) + sum(VyRim.*FyRim1));
    RecipSpeeds(i) = U;

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
