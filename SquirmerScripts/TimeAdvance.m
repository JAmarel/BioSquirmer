function [Ux_history, Uy_history, W_history, theta_history, x_cm_history, y_cm_history, separation_history, dt_history, time_history, ...
    x_history, y_history, Matrix_history,x_cm_history_Recip, y_cm_history_Recip,separation_history_Recip,speed_history_Recip, W_history_Recip, ...
    x_history_Recip, y_history_Recip]...
    = TimeAdvance(T, dt_o, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, R, a, Scale, B1, B2)
%Calls Solve_U at each increment to simulate time.

% T = Total simulation time
% dt_o = Time step increment
% xcoord, ycoord = Initial beast (blob) coordinates in enc frame
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

r_cm_history(1) = sqrt(xcoord(1)^2 + ycoord(1)^2);

% No data for these first entries.
%We start with this initial condition at i = 1, solve for parameters, and
%save them at i + 1.
Ux_history(1,:) = 0;
Uy_history(1,:) = 0;
W_history(1,:) = 0;
fx_history(1,:) = 0;
fy_history(1,:) = 0;
time_history(1,:) = 0;
separation_history(1,:) = R - sqrt(xcoord(1)^2 + ycoord(1)^2) - a;
separation_history(1,:) = R - abs(xcoord(1)); %Change for wall scenario



%% Lorentz Variables 
%This calls the Reciprocal theorem routine. This and the box labeled
%Lorentz below are separate from the Boundary element boxes. Should function
%singly.

x_cm_history_Recip = zeros([Steps+1,1]);
y_cm_history_Recip = zeros([Steps+1,1]);
r_cm_history_Recip = zeros([Steps+1,1]);
theta_history_Recip = zeros([Steps+1,1]);
x_history_Recip = zeros([Steps+1,length(xcoord)]);
y_history_Recip = zeros([Steps+1,length(xcoord)]);

separation_history_Recip = zeros([Steps+1,1]);
Ux_history_Recip = zeros([Steps+1,1]);
Uy_history_Recip = zeros([Steps+1,1]);
W_history_Recip = zeros([Steps+1,1]);
speed_history_Recip = zeros([Steps+1,1]);

theta_history_Recip(1) = theta_o;

VxRim_Recip = VxRim;
VyRim_Recip = VyRim;

xcoord_Recip = xcoord;
ycoord_Recip = ycoord;

theta_Recip = theta_o;

x_cm_history_Recip(1) = xcoord(1);
y_cm_history_Recip(1) = ycoord(1);
x_history_Recip(1,:) = xcoord;
y_history_Recip(1,:) = ycoord;
r_cm_history_Recip(1) = sqrt(xcoord(1)^2 + ycoord(1)^2);

separation_history_Recip(1) = R - sqrt(xcoord(1)^2 + ycoord(1)^2) - a;
Ux_history_Recip(1,:) = 0;
Uy_history_Recip(1,:) = 0;
W_history_Recip(1,:) = 0;
speed_history_Recip(1,:) = 0;


for i = 1:Steps
    
    %% Lorentz
    
    
%     [Ux_Recip, Uy_Recip, W_Recip] = Lorentz_solve_U(xcoord_Recip, ycoord_Recip, x_Enc, y_Enc, epsilon, VxRim_Recip, VyRim_Recip, NRim, theta_Recip);
%     
% %%% Scaling timesteps based on distance from enclosure.
% %%% Smaller steps when near the edge
%     r_cm_history_Recip(i) = sqrt(x_cm_history_Recip(i)^2 + y_cm_history_Recip(i)^2); %Beast CM
% %     if R - (r_cm_history_Recip(i) + a) < .5*a
% %         dt_Recip = dt_o/50;
% %     elseif R - (r_cm_history_Recip(i) + a) < 1.5*a
% %         dt_Recip = dt_o/10;
% %     elseif R - (r_cm_history_Recip(i) + a) < 2.5*a
% %         dt_Recip = dt_o/5;
% %     else
% %         dt_Recip = dt_o;
% %     end  
%     
% % Radial Velocity Check
% % V_r_Recip = (1/r_cm_history_Recip(i))*(Ux_Recip*x_cm_history_Recip(i) + Uy_Recip*y_cm_history_Recip(i));
% % %If beast would walk over the edge, change dt so that it is only possible
% % %to travel 1/4 of the distance remaining from the edge.
% %         if r_cm_history_Recip(i)+ a + V_r_Recip*dt_Recip > R
% %             dt_Recip = abs((1/(4*V_r_Recip))*(R - r_cm_history_Recip(i) - a));
% %         else
% %             dt_Recip = dt_Recip; %#ok<ASGSL>
% %         end
% %     
%     %%%Beast rotation
%     theta_Recip = theta_history_Recip(i) + W_Recip*dt_Recip;
% 
%     %Rotate for new velocities
%     [VxRim_Recip, VyRim_Recip] = UpdatedPrescribeWave(NRim, B1, B2, theta_Recip);
% 
% 
%     %Shift to beast frame (Translate origin to beast center)
%     xcoord_Recip = xcoord_Recip - x_cm_history_Recip(i);
%     ycoord_Recip = ycoord_Recip - y_cm_history_Recip(i);
% 
%     %Rotate beast coordinates due to W
%     [xcoord_Recip, ycoord_Recip] = Rotate_Vector(xcoord_Recip, ycoord_Recip, W_Recip*dt_Recip);
% 
%     %Back to Lab Frame.
%     xcoord_Recip = xcoord_Recip + x_cm_history_Recip(i);
%     ycoord_Recip = ycoord_Recip + y_cm_history_Recip(i);
% 
%     %%%Translation
%     xcoord_Recip = xcoord_Recip + Ux_Recip*dt_Recip;
%     ycoord_Recip = ycoord_Recip + Uy_Recip*dt_Recip;
%     
%     x_cm_history_Recip(i+1) = xcoord_Recip(1);
%     y_cm_history_Recip(i+1) = ycoord_Recip(1);
%     x_history_Recip(i+1,:) = xcoord_Recip;
%     y_history_Recip(i+1,:) = ycoord_Recip;
%     r_cm_history_Recip(i+1) = sqrt(xcoord_Recip(1)^2 + ycoord_Recip(1)^2);
%     
% 
%     separation_history_Recip(i+1,:) = R - (r_cm_history_Recip(i) + a);
%     Ux_history_Recip(i+1,:) = Ux_Recip;
%     Uy_history_Recip(i+1,:) = Uy_Recip;
%     W_history_Recip(i+1,:) = W_Recip;
%     speed_history_Recip(i+1,:) = sqrt(Ux_Recip^2 + Uy_Recip^2);
    
    
    
    
    
    
    
    
    %% Boundary Element
    %%% Solve for point forces, swimming velocity, angular rotation
    [fx, fy, Ux, Uy, W, Matrix, ~] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim, Scale);
% Step Forward
dt = dt_o;
Matrix_history{i} = Matrix;

  
%%% Scaling timesteps based on distance from enclosure.
%%% Smaller steps when near the edge
    r_cm_history(i) = sqrt(x_cm_history(i)^2 + y_cm_history(i)^2); %Beast CM
    separation_history(i) = R - (abs(x_cm_history(i))); %+a to get the blob nearest the edge
%     if separation_history(i) < .5*a
%         dt = dt_o/100;
%     elseif separation_history(i) < 1.5*a
%         dt = dt_o/25;
%     elseif separation_history(i) < 2.5*a
%         dt = dt_o/5;
%     else
%         dt = dt_o;
%     end
%     

% modified version for the wall
    %separation_history(i) = R - (x_cm_history(i) + a); %+a to get the blob nearest the edge
    if separation_history(i) < 0.2*a
        dt = dt_o/10;
    elseif separation_history(i) < 0.5*a
        dt = dt_o/5;
    elseif separation_history(i) < 1*a
        dt = dt_o/2;
    else
        dt = dt_o;
    end
    




     PercentCompleted = 100*i/Steps
    
    


  
% % Radial Velocity Check for boundary element
% V_r = (1/r_cm_history(i))*(Ux*x_cm_history(i) + Uy*y_cm_history(i));
% %If beast would walk over the edge, change dt so that it is only possible
% %to travel 1/4 of the distance remaining from the edge.
%         if r_cm_history(i)+ a + V_r*dt_o > R
%             dt = abs((1/(4*V_r))*(R - r_cm_history(i) - a));
%         else
%             dt = dt; %#ok<ASGSL>
%         end


            %%% If velocity unexpectedly changes by too much, go back to the previous step and
    %%% walk slightly further
    %%% history arrays are saved at i+1, so index i corresponds to the
    %%% previous timestep.
       if i>2 && sqrt(Ux^2 + Uy^2) > 3*sqrt(Ux_history(i)^2 + Uy_history(i)^2)
           Ux = Ux_history(i);
           Uy = Uy_history(i);
           W = W_history(i);
           dt = dt_history(i)/10;
       end




% % X component Velocity Check for wall situation
%If beast would walk over the edge, change dt so that it is only possible
%to travel 1/10 of the distance remaining from the edge.
        if x_cm_history(i)+ a + Ux*dt > R
            dt = abs((1/(10*Ux))*(R - x_cm_history(i) - a));
        end
        
        
        


        
            %For tracking timestep changes
            dt_history(i+1) = dt;
            time_history(i+1) = dt + sum(time_history(i));

            %%%Beast rotation
            theta = theta_history(i) + W*dt;


              
            %Rotate for new velocities
            [VxRim, VyRim] = UpdatedPrescribeWave(NRim, B1, B2, theta);


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
