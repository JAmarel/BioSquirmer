tic

%Simulation
T = 10;
dt = 1;


%Discretization
a = 0.1;              %%% radius of the disk nondimensionalized by the Saffman length
s= 0.1 * a;          %%% radial spacing between neighboring blobs
epsilon = s/8;       %%% radius of blobs

%Enclosure
R = 10*a;   %%%Radius of enclosure
d = 1*s;    %%%Circumferential Enclosure Blob Spacing

%Initial Conditions
r_o = -8.5*a;          %%% Radial coordinate of beast cm from center of enclosure
phi_o = 0*pi/4;     %%%Angle coordinate of beast cm from center of enclosure
theta_o = pi/2;   %%% Beast intial orientation (head direction)

%Coordinates of beast blobs in beast frame.
%Beast frame has its head on its x axis at y = 0.
[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

%Rotate coordinates according to theta_o
[xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, theta_o);

Nblobs = sum(BlobsPerLayer); %%% Number of blobs in the beast
NRim = BlobsPerLayer(end);   %%% Number of blobs in the outermost beast layer

%Enclosure Blob Coordinates from the origin in the enclosure frame.
[x_Enc, y_Enc] = DiscretizeEnclosure(R,d); 


%Prescribe wave in the beast frame.
[VxRim, VyRim, B1, B2] = PrescribeWave(NRim);

%Now Rotate velocities according to theta_o into lab frame.
[VxRim, VyRim] = Rotate_Vector(VxRim, VyRim, theta_o);

%Translate beast CM to (r_o, phi_o) in enclosure frame
x_o = r_o*cos(phi_o); %%%Beast CM initial x position as seen in enclosure frame.
y_o = r_o*sin(phi_o); %%%Beast CM Initial y position

xcoord = xcoord + x_o;
ycoord = ycoord + y_o;

%Grab the head coordinate for plotting purposes
x_head = xcoord(end - NRim + 1);
y_head = ycoord(end - NRim + 1);

[Ux_history, Uy_history, W_history, x_history, y_history, theta_history, x_cm_history, y_cm_history, fx_history, fy_history, COND_history] = ...
    TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, r_o, phi_o);


toc

%Create some strings for plot detail
str_T = ['T = ',num2str(T)];
str_dt = ['dt = ',num2str(dt)];
str_a = ['a = ',num2str(a)];
str_s = ['s = ',num2str(s)];
str_eps = ['epsilon = ',num2str(s)];
str_R = ['R = ',num2str(R)];
str_d = ['d = ',num2str(d)];
str_r_o = ['r_o = ',num2str(r_o)];
str_phi_o = ['phi_o = ',num2str(phi_o)];
str_theta_o = ['theta_o = ',num2str(theta_o)];
str_B1 = ['B1 = ',num2str(B1)];
str_B2 = ['B2 = ',num2str(B2)];
str_COND_max = ['COND_m_a_x = ', num2str(max(COND_history))];
str_COND_min = ['COND_m_i_n = ', num2str(min(COND_history))];

%% Plot the cm trajectory
fig = figure(1);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.3 .1 .6 .8]);
axes(ax1);
descr = {'Parameters:';
    str_T;
    str_dt;
    str_a;
    str_s;
    str_eps;
    str_R;
    str_d;
    str_r_o;
    str_phi_o;
    str_theta_o;
    str_B1;
    str_B2;
    str_COND_max;
    str_COND_min};
text(.025,0.6,descr)
%back to plotting data
axes(ax2);
plot(x_cm_history(1), y_cm_history(1), 'Marker', 'o', 'MarkerSize', 2, ...
     'MarkerFaceColor', 'green'); %Begin at green
hold on
plot(x_cm_history(end), y_cm_history(end), 'Marker', 'o', 'MarkerSize', 2, ...
     'MarkerFaceColor', 'red'); %End at red
plot(x_cm_history(2:end-1), y_cm_history(2:end-1),'Marker', '.', 'MarkerSize', 2, ...
     'MarkerFaceColor', 'black');
daspect([1,1,1])
%plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'r.','LineWidth', 1)
plot(x_Enc, y_Enc, 'Marker', '.', 'MarkerSize', 1, ...
     'MarkerFaceColor', 'black');
%plot(x_head, y_head, 'ko', 'LineWidth', 1)
axis off
hold off

%% Plot vector field at the last time step   %Dimensions may be wrong in here. 
% fig = figure(2);
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% ax2 = axes('Position',[.3 .1 .6 .8]);
% axes(ax1);
% descr = {'Parameters:';
%     str_T;
%     str_dt;
%     str_a;
%     str_s;
%     str_eps;
%     str_R;
%     str_d;
%     str_r_o;
%     str_phi_o;
%     str_theta_o;
%     str_B1;
%     str_B2};
% text(.025,0.6,descr)
% axes(ax2);
% rectangle('Position',[-R, -R, 2*R, 2*R],...
%            'Curvature',[1,1],...
%            'LineWidth', 2, 'LineStyle', '-', 'EdgeColor', 'k') %This outlines the enclosure
%        
% rectangle('Position',[x_cm_history(end) - a, y_cm_history(end) - a, 2*a, 2*a],...
%            'Curvature',[1,1],...
%            'LineWidth', 2, 'LineStyle', '-', 'EdgeColor', 'r') %This outlines the beast
% daspect([1,1,1])
% hold on
% 
% x = [x_cm_history(end) - 3*a :  0.15 * a : x_cm_history(end) + 3*a]';  %%% make a column
% y =  y_cm_history(end) - 3*a :  0.15 * a : y_cm_history(end) + 3*a;    %%% make a row
% 
% X = repmat(x, [1,length(y)]);   %%% form a matrix 
% Y = repmat(y, [length(x), 1]);  %%% form a matrix
% 
% VX = VX_FIELD_DISK(fx_history(end,:), fy_history(end,:), [x_Enc x_history(end,:)], [y_Enc y_history(end,:)], epsilon,  x, y);
% VY = VY_FIELD_DISK(fx_history(end,:), fy_history(end,:), [x_Enc x_history(end,:)], [y_Enc y_history(end,:)], epsilon,  x, y);
% 
% figure(2)
% quiver(X, Y, VX, VY, 'b')
% 
% xlim([x_cm_history(end) - 4*a x_cm_history(end) + 4*a]);
% ylim([y_cm_history(end) - 4*a y_cm_history(end) + 4*a]);
% 
% hold off