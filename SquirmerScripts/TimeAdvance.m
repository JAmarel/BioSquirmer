function [Ux_history, Uy_history, W_history, x_history, y_history, theta_history, Angles, x_cm_history, y_cm_history] = TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, r_o, phi_o)
%Calls Solve_U at each increment to simulate time.

% T = Total simulation time
% dt = Time steps
% xcoord, ycoord = Initial beast (blob) coordinates from origin
% x_Enc, y_Enc = Enclosure coodinates
% theta_o = Initial beast orientation

Steps = floor(T/dt); %Total number of increments

%Initialize the empty arrays
x_history = zeros([Steps+1,length(xcoord)]);
y_history = zeros([Steps+1,length(xcoord)]);

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


for i = 1:Steps
    [fx, fy, Ux, Uy, W, Ux_Enc, Uy_Enc, Matrix, N] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim);

    %%%Rotation
    theta = theta_o + W*dt;
    
    %Prescribe wave again to account for slight rotation.
    [VxRim, VyRim, B1] = PrescribeWave_Orient(NRim,theta);
 
    %First Shift to Beast Frame
    xcoord = xcoord - x_cm_history(i);
    ycoord = ycoord - y_cm_history(i);
    
    Angles = zeros([1, length(xcoord)]);
    Angles(1) = 0;
    %Get Angles for each blob from arctan. [-pi/2,pi/2]
    %These are the blob angles off the beast's x axis
    for j=2:length(xcoord)
        if xcoord(j)>=0 && ycoord(j)>=0
            Angles(j) = atan(abs(ycoord(j)/xcoord(j)));
        elseif xcoord(j)<0 && ycoord(j)>=0
            Angles(j) = pi - atan(abs(ycoord(j)/xcoord(j)));
        elseif xcoord(j)<0 && ycoord(j)<0
            Angles(j) = pi + atan(abs(ycoord(j)/xcoord(j)));
        elseif xcoord(j)>=0 && ycoord(j)<0
            Angles(j) = 2*pi - atan(abs(ycoord(j)/xcoord(j)));
        end
    end
    
    %New angles after W*dt
    NewAngles = Angles + W*dt;

    %Find new beast blob coordinates after rotation
    %Right now this is putting all blobs of one shell at the same coords.
    r_shell = sqrt(xcoord.^2 + ycoord.^2);
    xcoord = r_shell.*cos(NewAngles); %x(i+1) = r_shell*cos(theta(i+1))
    ycoord = r_shell.*sin(NewAngles);
    
    %Back to Lab Frame centered at the enclosure.
    xcoord = xcoord + x_cm_history(i);
    ycoord = ycoord + y_cm_history(i);
    
    %%%Translation
    xcoord = xcoord + Ux*dt;
    ycoord = ycoord + Uy*dt;
    
    %Fill in history
    x_history(i+1,:) = xcoord;
    y_history(i+1,:) = ycoord;
    
    x_cm_history(i+1) = xcoord(1);
    y_cm_history(i+1) = ycoord(1);
    
    Ux_history(i+1) = Ux;
    Uy_history(i+1) = Uy;

    W_history(1+1) = W;
    theta_history(1+1) = theta;
end

end

