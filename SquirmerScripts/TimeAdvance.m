function [Ux_history, Uy_history, W_history, x_history, y_history, theta_history, x_cm_history, y_cm_history, fx_history, fy_history] = TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, r_o, phi_o)
%Calls Solve_U at each increment to simulate time.

% T = Total simulation time
% dt = Time steps
% xcoord, ycoord = Initial beast (blob) coordinates in enclosure frame
% x_Enc, y_Enc = Enclosure coodinates
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

%Initial Positions
x_history(1,:) = xcoord;
y_history(1,:) = ycoord;

%CM is always first entry
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
    [fx, fy, Ux, Uy, W, Ux_Enc, Uy_Enc, Matrix, N] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim);

%     fx = fx/(B1/2); %Nondimensionalizing.
%     fy = fy/(B1/2);
%     Ux = Ux/(B1/2); %Now U and V have no dimensions.
%     Uy = Uy/(B1/2); %f carries dimensions [1/4pieta*h]
%     VxRim = VxRim/(B1/2);
%     VyRim = VyRim/(B1/2);
    
%     W = W/(B1/2);

    %%%Rotation
    theta = theta_o + W*dt;
    
    %rotate prescribed wave again account for W*dt.
    VxRimNew = VxRim*cos(W*dt) - VyRim*sin(W*dt);
    VyRimNew = VxRim*sin(W*dt) + VyRim*cos(W*dt);
    VxRim = VxRimNew;
    VyRim = VyRimNew;
 
    %Begin Shift to Beast Frame
    
    %First translate cm to origin
    xcoord = xcoord - x_cm_history(i);
    ycoord = ycoord - y_cm_history(i);
    
    %rotate coordinates to account for slight rotation.
    xcoordNew = xcoord*cos(W*dt) - ycoord*sin(W*dt);
    ycoordNew = xcoord*sin(W*dt) + ycoord*cos(W*dt);
    xcoord = xcoordNew;
    ycoord = ycoordNew;
    
    %Back to Lab Frame centered at the enclosure.
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
end

end

