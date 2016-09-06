 tic

%Simulation
T = 100;
dt = 1;

%Discretization
a = 1;              %%% radius of the disk nondimensionalized by the Saffman length
s = 0.1 * a;          %%% radial spacing between neighboring blobs
epsilon = s/8;        %%% radius of blobs
Scale = 10/a;

%Enclosure
R = 10*a;    %%% Radius of enclosure

%Initial Conditions
r_o = .35*R;         %%% Radial coordinate of beast cm from center of enclosure
phi_o = -pi/4;       %%% Angle coordinate of beast cm from center of enclosure
theta_o = .4*pi/2;   %%% Beast intial orientation (head direction)

x_o = r_o*cos(phi_o); %%% Beast CM initial x position as seen in enclosure frame.
y_o = r_o*sin(phi_o); %%% Beast CM Initial y position



%Coordinates of beast blobs in beast frame.
[xcoord, ycoord, BlobsPerLayer] = Discretize_Single_Shell_Disk(a,s);

Nbeast = sum(BlobsPerLayer); %%% Number of blobs in the beast
NRim = BlobsPerLayer(end);   %%% Number of blobs in the outermost beast layer

%Rotate coordinates according to theta_o
[xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, theta_o);

%Prescribe wave in the beast frame.
[VxRim, VyRim, B1, B2] = PrescribeWave(NRim);

%Now Rotate velocities according to theta_o into lab frame.
[VxRim, VyRim] = Rotate_Vector(VxRim, VyRim, theta_o);

%Nondimensionalize from the start?
VxRim = VxRim/(B1/2);
VyRim = VyRim/(B1/2);

%Enclosure Blob Coordinates.
%Enclosure is a ring centered about the lab origin.
[x_Enc, y_Enc] = DiscretizeEnclosure(R,s); 

%Translate beast geometric center to the chosen initial position in enclosure frame.
xcoord = xcoord + x_o;
ycoord = ycoord + y_o;

%Grab the head coordinate for plotting
x_head = xcoord(end - NRim + 1);
y_head = ycoord(end - NRim + 1);

[Ux_history, Uy_history, W_history, theta_history, x_cm_history, y_cm_history, separation_history, dt_history]...
    = TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, R, a, Scale);

speed_history = (Uy_history.^2 + Ux_history.^2).^(1/2);
NEnc = length(y_Enc);

toc

%% Create some strings for plot detail
str_T = ['Total Time = ',num2str(T)];
str_dt = ['Base dt = ',num2str(dt)];
str_Time = 'Nondimen by B_1/(2*l_s)';
str_a = ['a = ',num2str(a)];
str_s = ['s = ',num2str(s)];
str_eps = ['epsilon = ',num2str(epsilon)];
str_R = ['R = ',num2str(R)];
str_r_o = ['r_o = ',num2str(r_o)];
str_phi_o = ['phi_o = ',num2str(phi_o)];
str_theta_o = ['theta_o = ',num2str(theta_o)];
str_B1 = ['B1 = ',num2str(B1)];
str_B2 = ['B2 = ',num2str(B2)];
%str_COND_max = ['COND_m_a_x = ', num2str(max(COND_history))];
%str_COND_min = ['COND_m_i_n = ', num2str(min(COND_history))];

descr = {'Parameters:';
    str_T;
    str_dt;
    str_Time;
    str_a;
    str_s;
    str_eps;
    str_R;
    str_r_o;
    str_phi_o;
    str_theta_o;
    str_B1;
    str_B2};

%% Plot the cm trajectory
fig = figure(1);
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.3 .1 .6 .8]);
axes(ax1);
text(.025,0.6,descr);
%back to data
axes(ax2);

scatter(x_cm_history(1), y_cm_history(1), '.', 'g'); %Begin at green
hold on
scatter(x_cm_history(end-1), y_cm_history(end-1), '.', 'r'); %End at red
scatter(x_cm_history(2:end-2), y_cm_history(2:end-2), '.', 'k');
daspect([1,1,1])
%plot(xcoord(Nbeast - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'r.','LineWidth', 1)
scatter(x_Enc, y_Enc, '.', 'b');
%plot(x_head, y_head, 'ko', 'LineWidth', 1)
axis off
hold off
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