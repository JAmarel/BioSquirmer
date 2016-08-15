function [Ux_history, Uy_history, W_history, x_history, y_history, theta_history, x_cm_history, y_cm_history, fx_history, fy_history, COND_history, dt_history, r_cm_history, separation_history, time_history]...
    = TimeAdvance(T, dt_o, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, B1, R, a)
%Calls Solve_U at each increment to simulate time.

% T = Total simulation time
% dt = Time step increment
% xcoord, ycoord = Initial beast (blob) coordinates in lab frame
% x_Enc, y_Enc = Enclosure coodinates in lab frame
% theta_o = Initial beast orientation

Steps = floor(T/dt_o); %Total number of increments

%Initialize the empty arrays
dt_history = zeros([Steps+1,1]);
time_history = zeros([Steps+1,1]);

x_history = zeros([Steps+1,length(xcoord)]);
y_history = zeros([Steps+1,length(xcoord)]);

fx_history = zeros([Steps+1,length([xcoord x_Enc])]);
fy_history = zeros([Steps+1,length([xcoord x_Enc])]);

x_cm_history = zeros([Steps+1,1]);
y_cm_history = zeros([Steps+1,1]);

r_cm_history = zeros([Steps+1,1]);
separation_history = zeros([Steps+1,1]); %Distance between enclosure and the closest blob

Ux_history = zeros([Steps+1,1]);
Uy_history = zeros([Steps+1,1]);

W_history = zeros([Steps+1,1]);
theta_history = zeros([Steps+1,1]);

COND_history = zeros([Steps,1]);

%Initial Positions
x_history(1,:) = xcoord;
y_history(1,:) = ycoord;

%beast center is always first entry
x_cm_history(1) = xcoord(1);
y_cm_history(1) = ycoord(1);


theta_history(1) = theta_o;

% No data for these first entries.
Ux_history(1,:) = 0;
Uy_history(1,:) = 0;
W_history(1,:) = 0;
fx_history(1,:) = 0;
fy_history(1,:) = 0;
time_history(1,:) = 0;

dt = dt_o;

for i = 1:Steps
    [fx, fy, Ux, Uy, W, Matrix, N] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim);

    % Scale timesteps based on how close the enclosure is
    r_cm_history(i) = sqrt(x_cm_history(i)^2 + y_cm_history(i)^2);
    separation_history(i) = R - (r_cm_history(i) + a);
    
%%% Scaling timesteps based on distance from enclosure
%     if separation_history(i) < .5*a
%         dt = dt_o/50;
%     elseif separation_history(i) < 1*a
%         dt = dt_o/10;
%     else
%         dt = dt_o;
%     end

% Calculate the radial velocity
    V_r = (1/r_cm_history(i))*(Ux*x_cm_history(i) + Uy*y_cm_history(i));

%%% Limit radial distance traveled
%     if separation_history(i) < .5*a
%         dt = separation_history(i)/(4*;
%     elseif separation_history(i) < 1*a
%         dt = dt_o/10;
%     else
%         dt = dt_o;
%     end


    dt_history(i) = dt;
    time_history(i+1) = dt + sum(time_history(i));
    
    %%%Beast rotation
    theta = theta_o + W*dt;
    
    %Prescribe wave again account for W*dt.
    [VxRim, VyRim] = Rotate_Vector(VxRim, VyRim, W*dt);
 
    %Shift to beast frame
    %Translate origin at beast center
    xcoord = xcoord - x_cm_history(i);
    ycoord = ycoord - y_cm_history(i);
    
    %Rotate beast coordinates due to W
    [xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, W*dt);

    %Back to Lab Frame.
    xcoord = xcoord + x_cm_history(i);
    ycoord = ycoord + y_cm_history(i);
    
    %%%Translation
    xcoord = xcoord + Ux*dt;
    ycoord = ycoord + Uy*dt;
    
    %Fill in resulting (future) data
    x_history(i+1,:) = xcoord;
    y_history(i+1,:) = ycoord;
    
    fx_history(i+1,:) = fx;
    fy_history(i+1,:) = fy;
    
    x_cm_history(i+1) = xcoord(1);
    y_cm_history(i+1) = ycoord(1);
    
    Ux_history(i+1) = Ux;
    Uy_history(i+1) = Uy;

    W_history(i+1) = W;
    theta_history(i+1) = theta;
end  
