function [Ux_history, Uy_history, W_history, x_history, y_history, theta_history, Angles] = TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, r_o, phi_o)
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

Ux_history = zeros([Steps+1,1]);
Uy_history = zeros([Steps+1,1]);

W_history = zeros([Steps+1,1]);
theta_history = zeros([Steps+1,1]);

%Initial Positions
x_history(1,:) = xcoord;
y_history(1,:) = ycoord;
theta_history(1,:) = theta_o;

% No initial velocities
Ux_history(1,:) = 0;
Uy_history(1,:) = 0;
W_history(1,:) = 0;


for i = 1:Steps
    [fx, fy, Ux, Uy, W, Ux_Enc, Uy_Enc, Matrix, N] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim);

    %%%Rotation
    theta = theta_o + W*dt;
 
    %First Shift to Beast Frame
    xcoord = xcoord - r_o*cos(phi_o);
    ycoord = ycoord - r_o*sin(phi_o);
    
    Angles = zeros([1, length(xcoord)]);
    Angles(1) = 0;
    %Get Angles for each blob from arctan. [-pi/2,pi/2]
    %These are the blob angles off the beast's x axis
    for i=2:length(xcoord)
        if xcoord(i)>=0 && ycoord(i)>=0
            Angles(i) = atan(abs(ycoord(i)/xcoord(i)));
        elseif xcoord(i)<0 && ycoord(i)>=0
            Angles(i) = pi - atan(abs(ycoord(i)/xcoord(i)));
        elseif xcoord(i)<0 && ycoord(i)<0
            Angles(i) = pi + atan(abs(ycoord(i)/xcoord(i)));
        elseif xcoord(i)>=0 && ycoord(i)<0
            Angles(i) = 2*pi - atan(abs(ycoord(i)/xcoord(i)));
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
    xcoord = xcoord + r_o*cos(phi_o);
    ycoord = ycoord + r_o*sin(phi_o);
    
    %%%Translation
    xcoord = xcoord + Ux*dt;
    ycoord = ycoord + Uy*dt;
    
    %Fill in history
    x_history(i+1,:) = xcoord;
    y_history(i+1,:) = ycoord;
    
    Ux_history(i+1,:) = Ux;
    Uy_history(i+1,:) = Uy;

    W_history(1+1,:) = W;
    theta_history(1+1,:) = theta;
end

end

