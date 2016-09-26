function [fx_beast, fy_beast, fx_encl, fy_encl, Ux, Uy, Omega] = ...
    solve_levineslets_beast_encl(x0, y0, xbeast, ybeast, xencl, yencl, epsilon, VXbeast, VYbeast)



%x0 = x coordinate of beast center
%   This is xcoord(1)

%xbeast = array containing coordinates of all beast blobs. x0 first
%   This is xcoord

%xencl = array containing enclosure coordinates 
%   This is x_Enc

%VXbeast = Velocity of all beast blobs. [0 VxRim]
%   This is [zeros([Nbeast-NRim, 1]); VxRim];
%Nbeast = length(xcoord)
%NRim = length(VxRim)










Nbeast = length(xbeast);
Nencl = length(xencl);

VXbeast = VXbeast(:); %%% make sure it is a column
VYbeast = VYbeast(:); %%% make sure it is a column
xbeast = xbeast(:);         
ybeast = ybeast(:); 
xencl = xencl(:);        
yencl = yencl(:); 


%%% pre-compute distances between the blobs tiling the squirmer and the
%%% enclosure

distance11   = zeros([Nbeast, Nbeast]);  %%% distances between blobs on the squirmer
chihat11_x   = zeros([Nbeast, Nbeast]);  %%% x-comp of vector connecting i-th and j-th blobs
chihat11_y   = zeros([Nbeast, Nbeast]);  %%% y-comp of vector connecting i-th and j-th blobs
alpha_par11  = zeros([Nbeast, Nbeast]);
alpha_perp11 = zeros([Nbeast, Nbeast]);


for i = 1 : Nbeast  %%% runs over blobs tiling the beast
    for j = 1 : Nbeast  %%% runs over blobs
        if ( i ~= j )
            distance11(i,j) = sqrt((xbeast(i) - xbeast(j))^2 + (ybeast(i) - ybeast(j))^2);
            chihat11_x(i,j) = (xbeast(i) - xbeast(j))/distance11(i,j);
            chihat11_y(i,j) = (ybeast(i) - ybeast(j))/distance11(i,j);
 
            alpha_par11(i,j)  = mob_par_approx(distance11(i, j));
            alpha_perp11(i,j) = mob_perp_approx(distance11(i, j));  
        else
            alpha_par11(i,j) =  mob_par_approx(epsilon);
            alpha_perp11(i,j) = mob_perp_approx(epsilon);           
        end
    end
end

distance22   = zeros([Nencl, Nencl]);  %distances between blobs on enclosure
chihat22_x   = zeros([Nencl, Nencl]);  %%% x-comp of vector connecting i-th and j-th blobs
chihat22_y   = zeros([Nencl, Nencl]);  %%% y-comp of vector connecting i-th and j-th blobs
alpha_par22  = zeros([Nencl, Nencl]);
alpha_perp22 = zeros([Nencl, Nencl]);

for i = 1:Nencl
    for j=1:Nencl
         if ( i ~= j )
             distance22(i,j) = sqrt((xencl(i) - xencl(j))^2 + (yencl(i) - yencl(j))^2);
             chihat22_x(i,j) = (xencl(i) - xencl(j))/distance22(i,j);
             chihat22_y(i,j) = (yencl(i) - yencl(j))/distance22(i,j);
             alpha_par22(i,j)  = mob_par_approx(distance22(i, j));
             alpha_perp22(i,j) = mob_perp_approx(distance22(i, j));
         else
             alpha_par22(i,j) =  mob_par_approx(epsilon);
             alpha_perp22(i,j) = mob_perp_approx(epsilon);
         end
    end
end

distance12   = zeros([Nbeast, Nencl]);  %distances between blobs on squirmer and enclosure
chihat12_x   = zeros([Nbeast, Nencl]);  %%% x-comp of vector connecting i-th and j-th blobs
chihat12_y   = zeros([Nbeast, Nencl]);  %%% y-comp of vector connecting i-th and j-th blobs
alpha_par12  = zeros([Nbeast, Nencl]);
alpha_perp12 = zeros([Nbeast, Nencl]);
 
  for i = 1:Nbeast
     for j = 1:Nencl
         
         distance12(i,j) = sqrt((xbeast(i) - xencl(j))^2 + (ybeast(i) - yencl(j))^2);
         chihat12_x(i,j) = (xbeast(i) - xencl(j))/distance12(i,j);
         chihat12_y(i,j) = (ybeast(i) - yencl(j))/distance12(i,j);
         alpha_par12(i, j)  = mob_par_approx(distance12(i, j));
         alpha_perp12(i, j) = mob_perp_approx(distance12(i, j));       
                                 
     end
  end



%%% x-component of the vector field at point i of the squirmer due to
%%% x-Levineslet at point j of the squirmer


M11 = zeros([Nbeast, Nbeast]); 
  for i = 1:Nbeast
      for j = 1:Nbeast
          if ( i ~= j )
                M11(i, j) = chihat11_x(i, j)^2 * alpha_par11(i, j) + chihat11_y(i, j)^2 * alpha_perp11(i, j);                
          else
                M11(i, j) = alpha_par11(i, j);
          end
      end
  end
  

M21 = zeros([Nbeast, Nbeast]);  %%% y-component of the vector filed due to x-levinslet

  for i = 1:Nbeast
      for j = 1:Nbeast
          if (i ~= j )
                 M21(i,j) = chihat11_x(i,j) * chihat11_y(i,j) * ( alpha_par11(i, j) - alpha_perp11(i, j) );     
          else
                 M21(i,j) = 0;
          end
      end
  end

 
 
 M12 = M21;
  
 M22 = zeros([Nbeast, Nbeast]); %%% y-component of the velocity field due to y-levinslets
 for i = 1:Nbeast
      for j = 1:Nbeast
          if ( i ~= j ) 
                M22(i, j) = chihat11_y(i,j)^2 * alpha_par11(i, j) + ...
                            chihat11_x(i,j)^2 * alpha_perp11(i, j);                                     
             
          else
                M22(i, j) = alpha_perp11(i, j);
          end
      end
 end
 

%%% Consider contributions to the velocity field from levineslets on the
%%% enclosure


A11 = zeros([Nencl, Nencl]); %%% x-component of the vector field due to x-Levineslet
  for i = 1:Nencl
      for j = 1:Nencl
          if ( i ~= j )
                A11(i, j) = chihat22_x(i, j)^2 * alpha_par22(i, j) + chihat22_y(i, j)^2 * alpha_perp22(i, j);                  
          else
                A11(i, j) = alpha_par22(i, j);
          end
      end
  end
  
A21 = zeros([Nencl, Nencl]);  %%% y-component of the vector filed due to x-levineslet

  for i = 1:Nencl
      for j = 1:Nencl
          if (i ~= j )
                 A21(i,j) = chihat22_x(i,j) * chihat22_y(i,j) * ( alpha_par22(i, j) - alpha_perp22(i, j) );       
          else
                 A21(i,j) = 0;
          end
      end
  end

 A12 = A21;
  
 A22 = zeros([Nencl, Nencl]); %%% y-component of the velocity field due to y-levineslets
 for i = 1:Nencl
      for j = 1:Nencl
          if ( i ~= j ) 
                A22(i, j) = chihat22_y(i,j)^2 * alpha_par22(i, j) + ...
                            chihat22_x(i,j)^2 * alpha_perp22(i, j);                                   
             
          else
                A22(i, j) = alpha_perp22(i, j);
          end
      end
 end

%%%% Compute contributions to the velocity field at a point i on the squirmer due to
%%%% levineslets on the enclosure
  
  K11 = zeros([Nbeast, Nencl]); %%% x-component of the vector field on squirmer due to x-Levineslet on inlcusion enclosure
  
  for i = 1:Nbeast
      for j = 1:Nencl
          K11(i, j) = chihat12_x(i, j)^2 * alpha_par12(i, j) + chihat12_y(i, j)^2 * alpha_perp12(i, j);                
      end
  end
  
  K21 = zeros([Nbeast, Nencl]);  %%% y-component of the vector filed on squirmer due to x-Levineslet on inlcusion enclosure

  for i = 1:Nbeast
      for j = 1:Nencl
          K21(i,j) = chihat12_x(i,j) * chihat12_y(i,j) * ( alpha_par12(i, j) - alpha_perp12(i, j) );    
      end
  end

 K12 = K21;
  
 K22 = zeros([Nbeast, Nencl]); %%% y-component of the velocity field due to y-levinslets
 for i = 1:Nbeast
      for j = 1:Nencl
                K22(i, j) = chihat12_y(i,j)^2 * alpha_par12(i, j) + ...
                            chihat12_x(i,j)^2 * alpha_perp12(i, j);                                     
      end
 end
 
  
%%%% Consider now contributions to the velocity field on the enclosure from levineslets on
%%%% the squirmer

%%% Transpose matrices K11, K21, K22

    B11 = K11';  %%% B11 is Nencl x Nbeast matrix
    B21 = K21';  %%% transpose
    B12 = B21;
    B22 = K22';
  
  

 
 %%% Now we need to form one big matrix of the system of equations
 
 Matrix = [M11 M12 K11 K12 (-1)*ones([Nbeast,1])  zeros([Nbeast,1]) (-1)*(ybeast - y0); ...
           M21 M22 K21 K22  zeros([Nbeast, 1]) (-1)*ones([Nbeast,1]) (xbeast - x0)     ; ...
           B11 B12 A11 A12  zeros([Nencl,1])  zeros([Nencl,1]) zeros([Nencl,1]) ; ...
           B21 B22 A21 A22  zeros([Nencl,1])  zeros([Nencl,1]) zeros([Nencl,1]) ; ...                               
           ones([1, Nbeast]) zeros([1, Nbeast]) zeros([1,2*Nencl]) 0 0 0; ...               
           zeros([1, Nbeast])  ones([1, Nbeast])  zeros([1,2*Nencl])   0  0  0;...
            (-1)*(ybeast' - y0)  (xbeast' - x0) zeros([1,2*Nencl]) 0  0  0]; 
        
  c_coefs = [VXbeast; VYbeast; zeros([2*Nencl, 1]); 0; 0; 0];       
 

%%%Tests to debug 
%   Matrix = [M11 M12 K11 K12 (-1)*ones([Nbeast,1])  zeros([Nbeast,1]) ; ...
%            M21 M22 K21 K22  zeros([Nbeast, 1]) (-1)*ones([Nbeast,1])  ; ...
%            B11 B12 A11 A12  zeros([Nencl,1])  zeros([Nencl,1])  ; ...
%            B21 B22 A21 A22  zeros([Nencl,1])  zeros([Nencl,1]) ; ...                               
%            ones([1, Nbeast]) zeros([1, Nbeast]) zeros([1,2*Nencl]) 0 0 ; ...               
%            zeros([1, Nbeast])  ones([1, Nbeast])  zeros([1,2*Nencl])   0  0  ];
                  
        

 
%   Matrix = [ M11 M12 K11 K12 ; ...
%              M21 M22 K21 K22 ; ...
%              B11 B12 A11 A12  ; ...
%              B21 B22 A21 A22 ];                              
%                        
%           
%  c_coefs = [ones([Nbeast, 1]); zeros([Nbeast, 1]); zeros([2*Nencl, 1])];
%%%%%%%%%%%%%%

 variables = Matrix\c_coefs;
 
 fx_beast = variables(1:Nbeast);
 fy_beast = variables(Nbeast+1:2*Nbeast);
 fx_encl = variables(2*Nbeast + 1: 2*Nbeast + Nencl);
 fy_encl = variables(2*Nbeast + Nencl + 1: 2*Nbeast + 2*Nencl);
 Ux = variables(2*Nbeast + 2*Nencl + 1);
 Uy = variables(2*Nbeast + 2*Nencl + 2);
 Omega = variables(2*Nbeast + 2*Nencl + 3);

end