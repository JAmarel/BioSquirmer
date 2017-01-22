function [Ux, Uy, W] = ...
    Lorentz_solve_U(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim, theta_o)

%Forms the system of equations necessary to solve for the beast velocities
%in the Lorentz method. Each equation comes from calling
%enclosure_solve_U_disk once.



        
       

        %% Inactive Recip Thm
        
        % Inactive

        %Velocity of inactive disk
        U_hat = 1;
        W_hat = 1;
        
        
        Ux_hat = U_hat*cos(theta_o); 
        Uy_hat = U_hat*sin(theta_o);
        
        % Hat is conjugate scenario
        % Should I use the equations below, or FxNet_hat_1?
        % The value for Fx_hat is the resulting net force on the body necessary for our
        % choice of Ux_hat and Uy_hat according to the predictions of our
        % HPW mobility calculation.
        
        %Fx_hat = Ux_hat/HPW_mobility(a);
        %Fy_hat = Uy_hat/HPW_mobility(a);
        
        %The alternate choice is a similar object, the resulting net force
        %on the body necessary for our choice of Ux_hat and Uy_hat, but
        %this time according to the results of the Lorentz Reciprocal
        %theorem solver. I chose this option for now, because as we have to
        %use the rim forces from the lorentz solver, it seems best for
        %consistency to also use the total force of the same origin. The
        %values are both reasonable and similar.
        
        
        % Recip theorum no longer allows us to solve for Ux, Uy, and W separately.
        
        % First, solve for Ux as predicted by recip thm
        [fx_hat_1, fy_hat_1] = enclosure_solve_U_disk_inactive(xcoord, ycoord, x_Enc, y_Enc, epsilon, NRim, Ux_hat, 0, 0);

        FxRim_hat_1 = fx_hat_1(end-NRim+1:end);
        FyRim_hat_1 = fy_hat_1(end-NRim+1:end);
        
        
        %These form the first row of the 3x3 system
        Torque_hat_1 = dot(xcoord,fy_hat_1.') - dot(ycoord,fx_hat_1.');
        FxNet_hat_1 = sum(fx_hat_1); %%% x-component of net rim force on squirmer
        FyNet_hat_1 = sum(fy_hat_1); %%% y-component of net rim force on squirmer
        
        
       %Right hand side of Lorentz relationship. Answer column
        Sum_1 = -(sum(VxRim.*FxRim_hat_1) + sum(VyRim.*FyRim_hat_1));
        
        
        

        
        % Second, solve for Uy from recip thm
        [fx_hat_2, fy_hat_2] = enclosure_solve_U_disk_inactive(xcoord, ycoord, x_Enc, y_Enc, epsilon, NRim, 0, Uy_hat, 0);
        
        FxRim_hat_2 = fx_hat_2(end-NRim+1:end);
        FyRim_hat_2 = fy_hat_2(end-NRim+1:end);
        
        Torque_hat_2 = dot(xcoord,fy_hat_2.') - dot(ycoord,fx_hat_2.');
        FxNet_hat_2 = sum(fx_hat_2); %%% x-component of net rim force on squirmer
        FyNet_hat_2 = sum(fy_hat_2); %%% y-component of net rim force on squirmer
        
        Sum_2 = -(sum(VxRim.*FxRim_hat_2) + sum(VyRim.*FyRim_hat_2));
        
        % Third W case
        
        [fx_hat_3, fy_hat_3] = enclosure_solve_U_disk_inactive(xcoord, ycoord, x_Enc, y_Enc, epsilon, NRim, 0, 0, W_hat);
        
        FxRim_hat_3 = fx_hat_3(end-NRim+1:end);
        FyRim_hat_3 = fy_hat_3(end-NRim+1:end);
        
        Torque_hat_3 = dot(xcoord,fy_hat_3.') - dot(ycoord,fx_hat_3.');
        FxNet_hat_3 = sum(fx_hat_3); %%% x-component of net rim force on squirmer
        FyNet_hat_3 = sum(fy_hat_3); %%% y-component of net rim force on squirmer
        
        Sum_3 = -(sum(VxRim.*FxRim_hat_3) + sum(VyRim.*FyRim_hat_3));
        
        
        %Invert Matrix
        
        Matrix = [Torque_hat_1 FxNet_hat_1 FyNet_hat_1 ;...
          Torque_hat_2 FxNet_hat_2 FyNet_hat_2 ;... 
          Torque_hat_3 FxNet_hat_3 FyNet_hat_3];
      
       c_coefs = [Sum_1; Sum_2; Sum_3];
       
       %When the beast is facing 0,pi/2,pi,3pi/2, Sum_1 or 2 or 3 tend to 0
       %and make the equations noninvertible
       
       variables = Matrix\c_coefs;
       
       W = variables(1);
       Ux = variables(2);
       Uy = variables(3);
