function [Ux_history, Uy_history, W_history, x_history, y_history, theta_history, x_cm_history, y_cm_history, fx_history, fy_history, COND_history]...
    = TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, B1)
%Calls Solve_U at each increment to simulate time.

% T = Total simulation time
% dt = Time step increment
% xcoord, ycoord = Initial beast (blob) coordinates in lab frame
% x_Enc, y_Enc = Enclosure coodinates in lab frame
% theta_o = Initial beast orientation

Steps = floor(T/dt); %Total number of increments

%Initialize the empty arrays
x_history = zeros([Steps+1,length(xcoord)]);
y_history = zeros([Steps+1,length(xcoord)]);

fx_history = zeros([Steps+1,length([xcoord x_Enc])]);
fy_history = zeros([Steps+1,length([xcoord x_Enc])]);

x_cm_history = zeros([Steps+1,1]);
y_cm_history = zeros([Steps+1,1]);

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


for i = 1:Steps
    [fx, fy, Ux, Uy, W, Matrix, N] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim);

    %%%Beast rotation
    theta = theta_o + W*dt;
    
    %Prescribe wave again account for W*dt.
    [VxRim, VyRim] = Rotate_Vector(VxRim, VyRim, W*dt);
 
    %Shift to beast frame
    %Translate origin at beast center
    xcoord = xcoord - x_cm_history(i);
    ycoord = ycoord - y_cm_history(i);
    
    %Rotate beast coordinates to account for any rotation.
    [xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, W*dt);

    %Back to Lab Frame.
    xcoord = xcoord + x_cm_history(i);
    ycoord = ycoord + y_cm_history(i);
    
    %%%Translation
    xcoord = xcoord + Ux*dt;
    ycoord = ycoord + Uy*dt;
    
    %Fill in history
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
    
    COND_history(i) = cond(Matrix);
end

end

