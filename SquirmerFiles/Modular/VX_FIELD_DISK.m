function V = VX_FIELD_DISK(fx, fy, xcoord, ycoord, epsilon,  x, y)

%%% This routine computes x-component of the fluid velocity at positions
%%% xcoord and ycoord are arrays that contain coordinates of the blobs
%%% spread out on the disk
%%% fx and fy are x- and y- components of Levineslets distributed on the
%%% disk
%%% x and y are arrays that span the area over which we'd like to evaluate
%%% the velocity field; a pair (x(i), y(j)) gives coordinates of a point in x-y plane 

tiny = epsilon;  %%% watch for points very close to the island's circumference
L = length(xcoord); %%% total number of blobs

%%% pre-compute distances
distanceX = zeros([length(x), L]);
distanceY = zeros([length(y), L]);
distance = zeros([length(x), length(y), L]);

alpha_par = zeros([length(x), length(y), L]);
alpha_perp = zeros([length(x), length(y), L]);


%Euler = - psi(1);  
for i=1:length(x)
    for j=1:length(y)
        for k=1:L     %%% runs over blobs
           distanceX(i, k) = x(i) - xcoord(k);
           distanceY(j, k) = y(j) - ycoord(k);
           %%% compute distance between the k-th blob and a point in the
           %%% plane with coordinates given by x(i) and y(j)
           distance(i, j, k) = sqrt(distanceX(i, k)^2 + distanceY(j, k)^2);
       
           if (abs(distanceX(i, k)) > tiny || abs(distanceY(j,k)) > tiny )
           %%%if (distance(i, j, k) > tiny)  
              alpha_par(i, j, k) =  mob_par_approx(distance(i, j, k));
              alpha_perp(i, j, k) = mob_perp_approx(distance(i, j, k));
           else 
              alpha_par(i, j, k) = mob_par_approx(tiny);
              alpha_perp(i, j, k) =  mob_perp_approx(tiny);
           end
                   
        end
    end
end


Vx = zeros([length(x), length(y), L ]); %%% Vx-field due to x-levinslet at k-th blob
    
 for i = 1:length(x)
     for j = 1:length(y)
         for k = 1:L
             %%%if ( distance(i, j, k) > tiny ) 
             if (abs(distanceX(i, k)) > tiny || abs(distanceY(j,k)) > tiny ) 
                 Vx(i,j,k) = (1/distance(i, j, k))^2 * ...
                     ( (distanceX(i, k))^2 * alpha_par(i, j, k) + ...
                       (distanceY(j, k))^2 * alpha_perp(i, j, k )  );
             else
                 Vx(i,j,k) = alpha_par(i, j, k);
             end
              
          end
     end
 end
 Vx;     
        


 Wx = zeros([length(x), length(y), L ]);
 for i = 1:length(x)
     for j = 1:length(y)
         for k = 1:L
             %%%if ( distance(i, j, k) > tiny ) 
             if (abs(distanceX(i, k)) > tiny || abs(distanceY(j,k)) > tiny ) 
                 Wx(i,j,k) = (1/distance(i, j, k))^2 * distanceX(i, k) * distanceY(j, k) *...
                     (alpha_par(i, j, k) - alpha_perp(i, j, k ) );
             else
                 Wx(i,j,k) = 0;
             end
          end
     end
 end

 Wx ;
 

 for i = 1:length(x)
     for j = 1:length(y)
         sum1=0;
         sum2=0;
         for k = 1:L
             sum1 = sum1 + fx(k)* Vx(i, j, k);
             sum2 = sum2 + fy(k)* Wx(i, j, k);                 
         end
             V(i,j) = sum1 + sum2;   %%% gives x-component of the vector field at a point (xi, yj) due to both x- and y-levineslets
     end   
 end

end