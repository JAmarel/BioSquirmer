%% Single squirmer with rotation
tic
a = 10;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.08 * a;          %%% spacing between neighboring blobs
epsilon = s/8;       %%% radius of the blob

[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

Nblobs = sum(BlobsPerLayer); %%% total number of blobs 

NR = length(BlobsPerLayer); %%% Number of radial layers
NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

[VxRim, VyRim, B1] = PrescribeWave(NRim);

[fx, fy, Ux, Uy,W, Matrix] = solve_U_disk_rot(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);

fx = fx/(B1/2); %Nondimensionalizing.
fy = fy/(B1/2);
Ux = Ux/(B1/2); %Now U,W, and V have no dimensions.
Uy = Uy/(B1/2); %f carries dimensions [1/4pi eta]
                %those units cancel in efficiency calculation.
W = W/(B1/2);
VxRim = VxRim/(B1/2);
VyRim = VyRim/(B1/2);


FxRim = fx(end-NRim+1:end);
FyRim = fy(end-NRim+1:end);

FxNet = sum(fx); %%% x-component of net force on squirmer
FyNet = sum(fy); %%% y-component of net force on squirmer

TorqueNet = dot(xcoord,fy.') - dot(ycoord,fx.'); % Should be 0.
     
speed = sqrt(Ux^2 + Uy^2);

efficiency = CalcEfficiency(FxRim, FyRim, VxRim, VyRim, a, speed);

eigenvalues = eig(Matrix);



%%Plot the position blobs by xcoord and ycoord
% figure(1)
% plot(xcoord, ycoord, 'o')
% daspect([1,1,1])
% hold on
% plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'ro', 'LineWidth', 3)
% axis off
% hold off

%% plot the forces in the outer layer of the disk
% AnglesRim = zeros([1, NRim]);
% for i = 1:NRim
%     AnglesRim(i) = (i-1) * 2 * pi/NRim;
% end
% 
% figure(2)
% plot(AnglesRim, fx(Nblobs - NRim + 1:end), 'ro', 'LineWidth', 3)
% hold on
% plot(AnglesRim, fy(Nblobs - NRim + 1:end), 'bo', 'LineWidth', 3)
% hold off

%% Plot vector field   
%  figure(3)
%  rectangle('Position',[-a, -a, 2*a, 2*a],...
%            'Curvature',[1,1],...
%            'LineWidth', 2, 'LineStyle', '-', 'EdgeColor', 'r')
%  daspect([1,1,1])
%  hold on
% 
% x = [-2 * a :  0.15 * a : 2 * a]';  %%% make a column
% y =  -2 * a :  0.15 * a : 2 * a;    %%% make a row
% 
% X = repmat(x, [1,length(y)]);   %%% form a matrix 
% Y = repmat(y, [length(x), 1]);  %%% form a matrix
% 
% VX = VX_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x, y);
% VY = VY_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x, y);
% 
% figure(3)
% quiver(X, Y, VX, VY, 'b')
% hold off
% 
% %%% Plot fluid velocity along phi=0 direction
% span_r = linspace(0, 3*a, 30);
%   
%  vx = zeros([1,length(span_r)]);
%  vy = zeros([1,length(span_r)]);
%  x=vx;
%  y=vx;
%  phi=0;  %%% choose direction along the direction of net motion
% for i=1:length(span_r) 
%          x(i) = span_r(i) * cos(phi);
%          y(i) = span_r(i) * sin(phi);
%          vx(i) = VX_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x(i), y(i));
%          vy(i) = VY_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x(i), y(i));    
% end
% 
% figure(2)
% plot(span_r/a, vx/(B1/2), 'ro', 'LineWidth', 3)
% hold on
% plot(span_r/a, vy, 'bo', 'LineWidth', 3)
% hold off
%% plot the forces in different inner layers
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
%     BlobsPerLayer(NLayer);
%     length(fx(start+1:finish));
% 
%     figure(j)
%     plot(AnglesLayer, fx(start+1 : finish), 'ro', 'LineWidth', 2)
%     hold all
%     plot(AnglesLayer, fy(start+1:finish), 'bo', 'LineWidth', 2)
%     titlestr = strcat({'Forces on Blobs at NR = '},{' '},{num2str(j)});
%     title(titlestr);
%     xlabel('Blob Coordinate Angle (Radians)')
%     ylabel('Blob Force')
%     legend('f_x','f_y')
%     hold off
% end
toc
