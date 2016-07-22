function [fx, fy, Ux, Uy, W, Matrix, N] = ...
    solve_U_enclosure(xcoord, ycoord, x_Enc, y_Enc, epsilon, VxRim, VyRim, NRim)

%%% Now including rotation and enclosure. M = 7x7
%%% Relaxing the force constraint on the enclosure and remove 
%%%the ability for the enclosure
%%% to gain nonzero velocity

%%% angle = an array that contains instantenous angular position of
%%% point-like forces on the circumference
%%% radius = an array that contain instantenous radial position of
%%% point-like forces

%%% We assign the vector field distribution on the boundary of the circle


N = length([xcoord x_Enc]);   %%% number of blobs everywhere
%%% make sure these are columns
VxRim=VxRim(:);             %%% x-component of velocity at the rim of the disk                
VyRim=VyRim(:);             %%% y-component of velocity at the rim of the disk
NRim;                       %%% number of blobs on the rim of the disk
NBlobs = length(xcoord); %%% Number of beast blobs
NEnc = length(x_Enc);    %%% Number of enclosure blobs

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
Matrix = [M11 M12 [zeros([NEnc,1]); -ones([NBlobs,1])] zeros([N,1]) [zeros([NEnc,1]); ycoord.'] [-ones([NEnc,1]); zeros([NBlobs,1])] zeros([N,1]); ...
          M21 M22 zeros([N,1]) [zeros([NEnc,1]); -ones([NBlobs,1])] [zeros([NEnc,1]); -xcoord.'] zeros([N,1]) [-ones([NEnc,1]); zeros([NBlobs,1])]; ...
          -[zeros([1,NEnc]) ones([1, NBlobs])] zeros([1, N]) 0 0 0 0 0; ...
          zeros([1, N]) -[zeros([1,NEnc]) ones([1, NBlobs])] 0 0 0 0 0; ...
          [zeros([1,NEnc]) ycoord] [zeros([1,NEnc]) -xcoord] 0 0 0 0 0];
       
 %scondition = cond(Matrix);
 
 c_coefs = [zeros([N-NRim,1]); VxRim; zeros([N-NRim,1]); VyRim; 0; 0; 0];
 % [VxEnclosure; VxInner; VxRim; VyEnclosure; VyInner; VyRim; FxBeast; FyBeast; TorqueBeast; Fx_Enc, Fy_Enc]

 
 %[variables] = Gauss_method(2 * N + 2, Matrix, c_coefs);
 
 variables = Matrix\c_coefs;
 
 fx = variables(1:N);
 fy = variables(N+1:2*N);
 Ux = variables(2*N + 1);
 Uy = variables(2*N + 2);
 W  = variables(2*N + 3);

end