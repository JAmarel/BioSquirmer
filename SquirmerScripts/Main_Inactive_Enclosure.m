% Computes the beast swimming velocity when confined to an enclosure as predicted by the Lorentz
% reciprocal theorem at one instant in time

tic

B1 = 1;
B2 = 0;

%Discretization
a = 10;              %%% radius of the disk nondimensionalized by the Saffman length
s = 0.1 * a;          %%% radial spacing between neighboring blobs
epsilon = s/8;        %%% radius of blobs

%Enclosure
R = 5*a;    %%% Radius of enclosure
d = 10*s;    %%% Circumferential Enclosure Blob Spacing

%Initial Conditions
r_o = 0.5*R;         %%% Radial coordinate of beast cm from center of enclosure
phi_o = .8*pi;       %%% Angle coordinate of beast cm from center of enclosure
theta_o = pi/4;   %%% Beast intial orientation (head direction)

x_o = r_o*cos(phi_o); %%% Beast CM initial x position as seen in enclosure frame.
y_o = r_o*sin(phi_o); %%% Beast CM Initial y position

%Coordinates of beast blobs in beast frame.
[xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

Nblobs = sum(BlobsPerLayer); %%% Number of blobs in the beast
NRim = BlobsPerLayer(end);   %%% Number of blobs in the outermost beast layer

%Rotate coordinates according to theta_o
[xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, theta_o);

%Prescribe wave in the beast frame.
[VxRim, VyRim] = UpdatedPrescribeWave(NRim, B1, B2,theta_o);

%Enclosure Blob Coordinates.
%Enclosure is a ring centered about the lab origin.
[x_Enc, y_Enc] = DiscretizeEnclosure(R,d); 

%Translate beast geometric center to the chosen initial position in enclosure frame.
xcoord = xcoord + x_o;
ycoord = ycoord + y_o;

%Is this solving as it should? More inside
[Ux, Uy] = solve_U_inactive_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim,theta_o,a);

toc

%% Create some strings for plot detail
% str_T = ['Total Time = ',num2str(sum(dt_history))];
% str_dt = ['Base dt = ',num2str(dt)];
% str_Time = 'Nondimensionalized by B_1/(2*l_s)';
% str_a = ['a = ',num2str(a)];
% str_s = ['s = ',num2str(s)];
% str_eps = ['epsilon = ',num2str(s)];
% str_R = ['R = ',num2str(R)];
% str_d = ['d = ',num2str(d)];
% str_r_o = ['r_o = ',num2str(r_o)];
% str_phi_o = ['phi_o = ',num2str(phi_o)];
% str_theta_o = ['theta_o = ',num2str(theta_o)];
% str_B1 = ['B1 = ',num2str(B1)];
% str_B2 = ['B2 = ',num2str(B2)];
% str_COND_max = ['COND_m_a_x = ', num2str(max(COND_history))];
% str_COND_min = ['COND_m_i_n = ', num2str(min(COND_history))];
% 
% descr = {'Parameters:';
%     str_T;
%     str_dt;
%     str_Time;
%     str_a;
%     str_s;
%     str_eps;
%     str_R;
%     str_d;
%     str_r_o;
%     str_phi_o;
%     str_theta_o;
%     str_B1;
%     str_B2;
%     str_COND_max;
%     str_COND_min};

%% Plot the cm trajectory
% fig = figure(1);
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% ax2 = axes('Position',[.3 .1 .6 .8]);
% axes(ax1);
% text(.025,0.6,descr);
% %back to data
% axes(ax2);
% plot(x_cm_history(1), y_cm_history(1), 'Marker', 'o', 'MarkerSize', 2, ...
%      'MarkerFaceColor', 'green'); %Begin at green
% hold on
% plot(x_cm_history(end-1), y_cm_history(end-1), 'Marker', 'o', 'MarkerSize', 2, ...
%      'MarkerFaceColor', 'red'); %End at red
% plot(x_cm_history(2:end-2), y_cm_history(2:end-2),'Marker', '.', 'MarkerSize', 2, ...
%      'MarkerFaceColor', 'black');
% daspect([1,1,1])
% %plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'r.','LineWidth', 1)
% plot(x_Enc, y_Enc, 'Marker', '.', 'MarkerSize', 1, ...
%      'MarkerFaceColor', 'black');
% %plot(x_head, y_head, 'ko', 'LineWidth', 1)
% axis off
% hold off
%% Plotting velocity vs separation
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
%     str_B2;
%     str_COND_max;
%     str_COND_min};
% text(.025,0.6,descr);
% %back to plotting data
% axes(ax2);
% scatter(separation_history(2:end-1), Ux_history(2:end-1));
% hold on
% scatter(separation_history(2:end-1), Uy_history(2:end-1));
% legend('Ux','Uy','Location','Best')
% title('Swimming Speed vs. Distance from Enclosure Wall','FontSize',16,'FontWeight','bold')
% xlabel('Nondimensional Separation Distance (a/l_s)')
% ylabel('Nondimensional Swimming Speed')
% hold off
%% Plot vector field at the last time step   %Dimensions may be wrong in here. 
% fig = figure(3);
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% ax2 = axes('Position',[.3 .1 .6 .8]);
% axes(ax1);
% text(.025,0.6,descr);
% %back to plotting data
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
%% Plotting separation vs simulation time
% fig = figure(4);
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% ax2 = axes('Position',[.3 .1 .6 .8]);
% axes(ax1);
% text(.025,0.6,descr);
% %back to plotting data
% axes(ax2);
% scatter(time_history(2:end-1), separation_history(2:end-1));
% title('Distance from Enclosure Wall vs Time','FontSize',16,'FontWeight','bold')
% xlabel('Nondimensional Time')
% ylabel('Nondimensional Separation')