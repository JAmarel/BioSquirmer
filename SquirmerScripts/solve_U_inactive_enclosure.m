function [Ux, Uy] = ...
    solve_U_inactive_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim, theta,a)


% This is unused? Incomplete?

%%% Computes the beast swimming velocity as predicted by the Lorentz reciprocal theorum.
%%% Ux_hat and Uy_hat are the prescribed velocity components of the inactive conjugate swimmer. 

U_hat = 1; % This can be any number, choose 1.

%Inactive Initial Direction
Ux_hat = U_hat*cos(theta);
Uy_hat = U_hat*sin(theta);


N = length([xcoord x_Enc]);   %%% number of blobs everywhere

%%% make sure these are columns
VxRim=VxRim(:);             %%% x-component of velocity at the rim of the disk                
VyRim=VyRim(:);             %%% y-component of velocity at the rim of the disk
NRim;                       %%% number of blobs on the rim of the disk
NBlobs = length(xcoord); %%% Number of beast blobs
NEnc = length(x_Enc);    %%% Number of enclosure blobs


% xvect plays the same role as xcoord did before the enclosure was added
% This adds the enclosure blobs to the system of equations
% Now M11 and others can be calculated the same as before
xvect = [x_Enc xcoord]; %%% Vectors containing all blob coordinates
yvect = [y_Enc ycoord];


%%% precompute distances between points i and j
distance   = zeros([N, N]);  %%% distance between i-th and j-th blob
alpha_par  = zeros([N, N]);
alpha_perp = zeros([N, N]);
chihat_x = zeros([N, N]);    %%% x-component of unit vector from i-th blob to j-th blob
chihat_y = zeros([N, N]);    %%% y-component of unit vector from i-th blob to j-th blob


%Euler = - psi(1); 
for i = 1 : N  %%% Pick a blob
    for j = 1 : N  %%% run over all other blobs
        if ( i ~= j )
            distance(i,j) = sqrt((xvect(i) - xvect(j))^2 + (yvect(i) - yvect(j))^2);
            chihat_x(i,j) = (xvect(i) - xvect(j))/distance(i,j);
            chihat_y(i,j) = (yvect(i) - yvect(j))/distance(i,j);

            alpha_par(i,j)  = mob_par_approx(distance(i, j));
            alpha_perp(i,j) = mob_perp_approx(distance(i, j));
        else
            
            alpha_par(i,j) =  mob_par_approx(epsilon);
            alpha_perp(i,j) = mob_perp_approx(epsilon);
        end
    end
end


%%% x-component of the vector field at blob i on the disk due to
%%% x-Levineslet at blob j on the circle


M11 = zeros([N, N]); 
  for i = 1:N
      for j = 1:N
          if ( i ~= j )
                M11(i, j) = chihat_x(i,j)^2 * alpha_par(i, j) + chihat_y(i,j)^2 * alpha_perp(i , j); 
          else
                M11(i, j) = alpha_par(i, j);
          end
      end
  end



M21 = zeros([N, N]);  %%% y-component of the vector field at point i due to x-levinslet at j

  for i = 1:N
      for j = 1:N
          if (i ~= j )

                 M21(i,j) = chihat_x(i,j) * chihat_y(i,j) * ( alpha_par(i, j) - alpha_perp(i, j) ); 
          else
                 M21(i,j) = 0;
          end
      end
  end


 M12 = M21;
  
 M22 = zeros([N, N]); %%% y-component of the velocity field due to y-levinslets
 for i = 1:N
      for j = 1:N
          if ( i ~= j ) 
                       
                M22(i, j) = chihat_y(i,j)^2 * alpha_par(i, j) + ...
                            chihat_x(i,j)^2 * alpha_perp(i, j);            
             
          else
                
                M22(i, j) = alpha_par(i, j);
          end
      end
 end
 

 %%% y-component of the velocity field at point i on the circumference due to a y-levinslet
 %%% in the center of the circle
 
 %%% Now we need to form one big matrix for the system of equations
 
%format short
Matrix = [M11 M12; ...
          M21 M22];
       
 %scondition = cond(Matrix);
 
 c_coefs = [zeros([NEnc,1]); Ux_hat*ones([NBlobs,1]); zeros([NEnc,1]); Uy_hat*ones([NBlobs,1])];
 % [[VxEnclosure; Ux_hat]; [VyEnclosure; Uy_hat]]

 
 variables = Matrix\c_coefs;
 
 fx_hat = variables(1:N);
 fy_hat = variables(N+1:2*N);
 
 fx_hat_rim = fx_hat(end-NRim+1:end);
 fy_hat_rim = fy_hat(end-NRim+1:end);
 
 
 % There is some confusion here on whether to use Fx_hat or F_hat_x as
 % calculated below.
 
 %Think this option is incorrect
 Fx_hat = sum(fx_hat);
 Fy_hat = sum(fy_hat);
 
 %Need to implement this option
 F_hat_x = Ux_hat/HPW_mobility(a);
 F_hat_y = Uy_hat/HPW_mobility(a);
 
 %How to correctly solve for Ux, Uy, and W after introducing the enclosure?
 %Can I still split them up? as in Main_Inactive?
 %Another system of equations?
 
 % Now use the reciprocal theorem to determine the active swimming speed
 
 Ux = (-1/Fx_hat)*(sum(VxRim.*fx_hat_rim) + sum(VyRim.*fy_hat_rim));
 Uy = (-1/Fy_hat)*(sum(VxRim.*fx_hat_rim) + sum(VyRim.*fy_hat_rim));
end