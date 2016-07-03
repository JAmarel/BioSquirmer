function MAIN
tic
a = 10;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.1 * a;          %%% spacing between neighboring blobs
epsilon = s/8;       %%% radius of the blob
NR = floor(a/s) + 1;  %%% number of layers along radial direction (+1 b/c of one blob that sits in the center)
B1 = 0.1;            %%% tangential velocity strength (streaming)
B2 = 0;   

%%% floor rounds the number to the nearest integer less or equal to that number

BlobsPerLayer = zeros([1, length(NR)]);

BlobsPerLayer(1) = 1; %%% one blob in the center of the circle
Nblobs = 1;
for i = 2: NR
    BlobsPerLayer(i) = round( 2 * pi * (i - 1) ); 
    Nblobs = Nblobs + BlobsPerLayer(i);
end

Nblobs; %%% total number of blobs 
NRim = BlobsPerLayer(NR)  %%% number of blobs in the outermost layer

%%% we need to find the coordinates of the blobs
xcoord = zeros([1, Nblobs]);
ycoord = xcoord;

xcoord(1)=0; %%% blob in the center of the circle
ycoord(1)=0; %%% blob in the center of the circle

index = 2;
for i=2:NR
    for j=1:BlobsPerLayer(i)
        delta_phi = 2*pi/BlobsPerLayer(i);
        xcoord(index) = (i-1) * s * cos((j-1) * delta_phi); 
        ycoord(index) = (i-1) * s * sin((j-1) * delta_phi);
        index = index + 1;
    end
end

figure(1)
plot(xcoord, ycoord, 'o')
daspect([1,1,1])
hold on
plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'ro', 'LineWidth', 3)
axis off
hold off



%%% prescribe tangential velocity on the rim of the disk (only)
VRimTheta = zeros([NRim, 1]); %%% tangential velocity at the rim of the disk
VxRim = zeros([NRim, 1]);     %%% x-compnent of tangential velocity at the rim
VyRim = zeros([NRim, 1]);     %%% y-compnent of tangential velocity at the rim

for i=1:NRim
    angle = (i-1) * 2 * pi/NRim;
    VRimTheta(i) = B1 * sin(angle) + B2 * sin(2 * angle);
    VxRim(i) = -VRimTheta(i) * sin(angle); 
    VyRim(i) = VRimTheta(i) * cos(angle);
end

       
%%% we use Levine and McKintosh response function
%%% v_\alpha = \alpha_{\alpha\beta} f_\beta
%%% in the code the response function is non-dimensionalized by multiplying
%%% it by 4\pi\eta h, where \eta is (3D) viscosity of the membrane and h is
%%% the thickness of the membrane: \tilde\alpha_\alpha\beta = alpha_\alpha\beta * 4\pi \eta h 
%%% Thus, in the code we actually use v_\alpha = \tilde\alpha_{\alpha\beta}
%%% * \tilde{f}_\beta, with \tilde{f}_\beta= \equiv f_\beta/ (4\pi\eta h),
%%% and \tilde{f}_\beta has units of m/s.


%%% Working in the lab frame we compute the forces on the squirmer as well
%%% as the translational velocity Ux and Uy 


[fx, fy, Ux, Uy] = ...
    solve_U_disk(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);
       
 FxNet = sum(fx) %%% x-component of net force on squirmer
 FyNet = sum(fy) %%% y-component of net force on squirmer

 Ux;
 Uy;
     
speed = Ux/(B1/2)   %%%% swimming velocity non-dimensionalized by B1/2

%% plot the forces in the outer layer of the disk
% AnglesRim = zeros([1, NRim]);
% for i = 1:NRim
%     AnglesRim(i) = (i-1) * 2 * pi/NRim;
% end
% 
% figure(4)
% plot(AnglesRim, fx(Nblobs - NRim + 1:end), 'ro', 'LineWidth', 3)
% hold on
% plot(AnglesRim, fy(Nblobs - NRim + 1:end), 'bo', 'LineWidth', 3)
% hold off

%% plot the vector field in different inner layers
% for j = 1:NR-1 %%% all inner layers
%     NLayer = j; %%% layer number for which we'd like to plot the forces
%     AnglesLayer = zeros([1, BlobsPerLayer(NLayer)]);
%     for i = 1:BlobsPerLayer(NLayer)
%         AnglesLayer(i) = (i-1) * 2 * pi/BlobsPerLayer(NLayer);
%     end
%     %%% 
%     start = 0;
%     for i = 1: NLayer - 1
%         start = start + BlobsPerLayer(i);    
%     end
%     finish = start + BlobsPerLayer(NLayer);
% 
%     BlobsPerLayer(NLayer)
%     length(fx(start+1:finish))
% 
%     figure(j)
%     plot(AnglesLayer, fx(start+1 : finish), 'ro', 'LineWidth', 3)
%     hold all
%     plot(AnglesLayer, fy(start+1:finish), 'bo', 'LineWidth', 3)
%     hold off
%     
% end

%% Plot vector field   
 figure(3)
 rectangle('Position',[-a, -a, 2*a, 2*a],...
           'Curvature',[1,1],...
           'LineWidth', 2, 'LineStyle', '-', 'EdgeColor', 'r')
 daspect([1,1,1])
 hold on

x = [-2 * a :  0.15 * a : 2 * a]';  %%% make a column
y =  -2 * a :  0.15 * a : 2 * a;    %%% make a row

X = repmat(x, [1,length(y)]);   %%% form a matrix 
Y = repmat(y, [length(x), 1]);  %%% form a matrix

VX = VX_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x, y);
VY = VY_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x, y);

figure(3)
quiver(X, Y, VX, VY, 'b')
hold off

%%% Plot fluid velocity along phi=0 direction
span_r = linspace(0, 3*a, 30);
  
 vx = zeros([1,length(span_r)]);
 vy = zeros([1,length(span_r)]);
 x=vx;
 y=vx;
 phi=0;  %%% choose direction along the direction of net motion
for i=1:length(span_r) 
         x(i) = span_r(i) * cos(phi);
         y(i) = span_r(i) * sin(phi);
         vx(i) = VX_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x(i), y(i));
         vy(i) = VY_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x(i), y(i));    
end

figure(2)
plot(span_r/a, vx/(B1/2), 'ro', 'LineWidth', 3)
hold on
plot(span_r/a, vy, 'bo', 'LineWidth', 3)
hold off
% toc
end