function [Ux_history, Uy_history, W_history, theta_history, x_cm_history, y_cm_history, separation_history, dt_history, time_history, ...
    x_history, y_history]...
    = TimeAdvance_Wall(T, dt_o, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, R, a, Scale)
%Calls Solve_U at each increment to simulate time.

% T = Total simulation time
% dt_o = Time step increment
% xcoord, ycoord = Initial beast (blob) coordinates in enc frame
% x_Enc, y_Enc = Wall coodinates in lab frame
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

separation_history = zeros([Steps+1,1]); %Distance between enclosure and the closest beast blob

Ux_history = zeros([Steps+1,1]);
Uy_history = zeros([Steps+1,1]);

W_history = zeros([Steps+1,1]);
theta_history = zeros([Steps+1,1]);

% For calculating the condition number at each step. Comp Intensive.
% COND_history = zeros([Steps,1]);

%Beast Initial Positions
x_history(1,:) = xcoord;
y_history(1,:) = ycoord;

%beast center is always first entry
x_cm_history(1) = xcoord(1);
y_cm_history(1) = ycoord(1);


theta_history(1) = theta_o;
separation_history(1,:) = R - abs(xcoord(1)); %Change for wall scenario

% No data for these first entries.
%We start with this initial condition at i = 1, solve for parameters, and
%save them at i + 1.
Ux_history(1,:) = 0;
Uy_history(1,:) = 0;
W_history(1,:) = 0;
fx_history(1,:) = 0;
fy_history(1,:) = 0;
time_history(1,:) = 0;


for i = 1:Steps 
    
    

    
    %% Boundary Element
    %%% Solve for point forces, swimming velocity, angular rotation
    %%% Enclosure solver can handle the wall just the same
    [fx, fy, Ux, Uy, W, ~, ~] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim, Scale);
    
    % Begin step forward with new velocities



  

    % Scale timesteps as the swimmer approaches the wall.
    separation_history(i) = R - x_cm_history(i); %Distance from beast center to wall
    if separation_history(i) < (1+0.2*a) %If the nearest edge of the squirmer is within 0.2 a of the wall.
        dt = dt_o/10;
    elseif separation_history(i) < (1+0.5*a)
        dt = dt_o/5;
    elseif separation_history(i) < (1+1*a)
        dt = dt_o/2;
    else
        dt = dt_o;
    end
    
    %If beast would walk over the edge, change dt so that it is only possible
    %to travel 1/4 of the distance remaining from the edge.
        if xcoord(1) + a + Ux*dt > R
            dt = abs((1/(4*Ux))*(separation_history(i) - a));
        else
            dt = dt; %#ok<ASGSL>
        end


     PercentCompleted = 100*i/Steps
    
        
            %For tracking timestep changes
            dt_history(i+1) = dt;
            time_history(i+1) = dt + sum(time_history(i));

            %%%Beast rotation
            theta = mod(theta_history(i) + W*dt, 2*pi); %When converting to Crowdy coordinate system, modulo addition is helpful.
              
            %Rotate rim velocities accordingly
            [VxRim, VyRim] = Rotate_Vector(VxRim, VyRim, W*dt);

            %Shift to beast frame (Translate origin to beast center)
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

            %Fill in future data
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
