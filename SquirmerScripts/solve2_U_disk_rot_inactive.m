function [fx, fy, Ux, Uy, W, Matrix] = ...
    solve2_U_disk_rot_inactive(xcoord, ycoord, epsilon, NRim, FxBeast, FyBeast)
%Inactive disk. Fnet !=0. Input Vrim = 0.
%%% Now including rotation. M = 5x5

%%% angle = an array that contains instantenous angular position of
%%% point-like forces on the circumference
%%% radius = an array that contain instantenous radial position of
%%% point-like forces

%%% We assign the vector field distribution on the boundary of the circle


N = length(xcoord);   %%% number of blobs in the disk
NRim;                       %%% number of blobs on the rim of the disk


VxRim = zeros([NRim,1]);
VyRim = zeros([NRim,1]);




%%% precompute distances between points i and j in the disk 
distance   = zeros([N, N]);  %%% distance between i-th and j-th blob
alpha_par  = zeros([N, N]);
alpha_perp = zeros([N, N]);
chihat_x = zeros([N, N]);    %%% x-component of unit vector from i-th blob to j-th blob
chihat_y = zeros([N, N]);    %%% y-component of unit vector from i-th blob to j-th blob


%Euler = - psi(1); 
for i = 1 : N  %%% runs over blobs, the last Nrim blobs are on the disks's rim
    for j = 1 : N  %%% runs over blobs
        if ( i ~= j )
            distance(i,j) = sqrt((xcoord(i) - xcoord(j))^2 + (ycoord(i) - ycoord(j))^2);
            chihat_x(i,j) = (xcoord(i) - xcoord(j))/distance(i,j);
            chihat_y(i,j) = (ycoord(i) - ycoord(j))/distance(i,j);

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
 
 %%% Now we need to form one big matrix of the system of equations
 
%format short       
Matrix = [M11 M12 -ones([N,1]) zeros([N,1]) ycoord.';...
          M21 M22 zeros([N,1]) -ones([N,1]) -xcoord.'; ...
          ones([1, N]) zeros([1, N]) 0 0 0;...
          zeros([1, N]) ones([1, N]) 0 0 0;...
          ycoord -xcoord 0 0 0];

       
 %scondition = cond(Matrix);
 

 c_coefs = [zeros([N-NRim,1]); VxRim; zeros([N-NRim,1]); VyRim; FxBeast; FyBeast; 0];
 % [VxInner;VxRim;VyInner;VyRim;FxBeast;FyBeast, Torque]
 
 %[variables] = Gauss_method(2 * N + 2, Matrix, c_coefs);
 
 variables = Matrix\c_coefs;
 
 fx = variables(1:N);
 fy = variables(N+1:2*N);
 Ux = variables(2*N + 1);
 Uy = variables(2*N + 2);
 W  = variables(2*N + 3);


end