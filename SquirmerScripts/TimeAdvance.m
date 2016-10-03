function [Ux_history, Uy_history, W_history, theta_history, x_cm_history, y_cm_history, separation_history, dt_history, time_history, x_history,Matrix_history]...
    = TimeAdvance(T, dt_o, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, R, a, Scale, B1, B2)
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


%%% Change to compare with Dr. K
% fx_history = zeros([Steps+1,length([xcoord])]);
% fy_history = zeros([Steps+1,length([xcoord])]);
%%% End Change


x_cm_history = zeros([Steps+1,1]);
y_cm_history = zeros([Steps+1,1]);

r_cm_history = zeros([Steps+1,1]);
separation_history = zeros([Steps+1,1]); %Distance between enclosure and the closest blob

Ux_history = zeros([Steps+1,1]);
Uy_history = zeros([Steps+1,1]);

W_history = zeros([Steps+1,1]);
theta_history = zeros([Steps+1,1]);

% Previously, needed to calculate the condition number at each step
% COND_history = zeros([Steps,1]);

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


for i = 1:Steps
    
    %%% Comment to compare
    [fx, fy, Ux, Uy, W, Matrix, ~] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim, Scale);
    %%% End
    



%%% For comparing with Dr. Kuriabova.
%%% Be sure to toggle the definitions for fx and fy found near the top.
% Nbeast = length(xcoord);
% VXbeast = [zeros([Nbeast-NRim, 1]); VxRim];
% VYbeast = [zeros([Nbeast-NRim, 1]); VyRim];
% x0 = xcoord(1);
% y0 = ycoord(1);
% xbeast = xcoord;
% ybeast = ycoord;
% xencl = x_Enc;
% yencl = y_Enc;
% 
% [fx_beast, fy_beast, fx_encl, fy_encl, Ux, Uy, Omega, Matrix] = ...
%     solve_levineslets_beast_encl(x0, y0, xbeast, ybeast, xencl, yencl, epsilon, VXbeast, VYbeast);
% fx = fx_beast;
% fy = fy_beast;
% W = Omega;
%%% End Change

Matrix_history{i} = Matrix;

  
%%% Scaling timesteps based on distance from enclosure
    r_cm_history(i) = sqrt(x_cm_history(i)^2 + y_cm_history(i)^2);
    separation_history(i) = R - (r_cm_history(i) + a);
    if separation_history(i) < .5*a
        dt = dt_o/100;
    elseif separation_history(i) < 1.5*a
        dt = dt_o/10;
    elseif separation_history(i) < 2.5*a
        dt = dt_o/2;
    else
        dt = dt_o;
    end
    
    PercentCompleted = 100*i/Steps
    
    
    %%% If velocity unexpectedly rises, make the next dt smaller
    %%% history arrays are saved at i+1, so index i corresponds to the
    %%% previous timestep.
%        if i>2 && sqrt(Ux^2 + Uy^2) > 3*sqrt(Ux_history(i)^2 + Uy_history(i)^2)
%            Ux = Ux_history(i);
%            Uy = Uy_history(i);
%            W = W_history(i);
%            dt = dt_history(i)/10;
%              ratio = sqrt(Ux_history(i)^2 + Uy_history(i)^2) / sqrt(Ux^2 + Uy^2);
%              dt = ratio*dt;
%        end

            
% Radial Velocity
V_r = (1/r_cm_history(i))*(Ux*x_cm_history(i) + Uy*y_cm_history(i));
%If I would increment over the edge, change dt so that it is only possible
%to travel 1/4 of the distance remaining from the edge.
        if r_cm_history(i)+ a + V_r*dt_o > R
            dt = (1/(4*V_r))*(R - r_cm_history(i) - a);
        else
            dt = dt_o;
        end
        
            dt_history(i+1) = dt;
            time_history(i+1) = dt + sum(time_history(i));

            %%%Beast rotation
            theta = theta_history(i) + W*dt;


              

            [VxRim, VyRim] = UpdatedPrescribeWave(NRim, B1, B2, theta);


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
