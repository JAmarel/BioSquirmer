% First inactive attempt. Zero prescribed wave velocity and nonzero
% prescribed body force. Choose an F and solve for U. This approach
% is dropped in favor of Pick a U and solve for F.

a = 10;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.08 * a;          %%% spacing between neighboring blobs
epsilon = s/8;       %%% radius of the blob

%Net body forces on beast
FxBeast = 8; 
FyBeast = 0;

[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

Nblobs = sum(BlobsPerLayer); %%% total number of blobs 

NR = length(BlobsPerLayer); %%% Number of radial layers
NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

[VxRim, VyRim, B1] = PrescribeWave(NRim);

%Is this the best way to nondimensionalize?
VxRim = VxRim/(B1/2);
VyRim = VyRim/(B1/2);

%% Inactive
[fx1, fy1, Ux1, Uy1, W1, Matrix1] = solve2_U_disk_rot_inactive(xcoord, ycoord, epsilon, NRim, FxBeast, FyBeast);

FxRim1 = fx1(end-NRim+1:end);
FyRim1 = fy1(end-NRim+1:end);

FxNet1 = sum(fx1); %%% x-component of net force on squirmer
FyNet1 = sum(fy1); %%% y-component of net force on squirmer

TorqueNet1 = dot(xcoord,fy1.') - dot(ycoord,fx1.');
     
speed1 = sqrt(Ux1^2 + Uy1^2);

%% Active
[fx2, fy2, Ux2, Uy2, W2, Matrix2] = solve_U_disk_rot(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);

FxRim2 = fx2(end-NRim+1:end);
FyRim2 = fy2(end-NRim+1:end);

FxNet2 = sum(fx2); %%% x-component of net force on squirmer
FyNet2 = sum(fy2); %%% y-component of net force on squirmer

TorqueNet2 = dot(xcoord,fy2.') - dot(ycoord,fx2.');
     
speed2 = sqrt(Ux2^2 + Uy2^2);
% 
%% Recip Thm

F_hat_x = Ux1/HPW_mobility(a);
F_hat_y = Uy1/HPW_mobility(a);

% F_hat_x = FxBeast;
% F_hat_y = FyBeast;

% These terms are from equation 12 in our notes
LHS = F_hat_x*Ux2 + F_hat_y*Uy2
RHS = -(sum(VxRim.*FxRim1) + sum(VyRim.*FyRim1))

Ratio = LHS/RHS % This ratio is 1 if we are in perfect agreement.