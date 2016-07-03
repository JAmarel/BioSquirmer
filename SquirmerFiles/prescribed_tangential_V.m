function prescribed_tangential_V

%%% semi-axes of an ellipse
%%% a - along x-axis and b-along y-axis

index = -2:0.25:5;
a = 10.^index;         %%% in ls units
B1 = 10;               %%% tangential velocity strength (streaming)
B2 = 0;                %%% tangential velocity strength (pushing)
N = 75;               %%% number of Levineslets on the boundary of squirmer
  

angles = linspace(0, 2*pi - 2*pi/N, N);
length(angles);
Del = 2*pi/N;           %%% Del is the angular segment
 
Ux = zeros([1,length(a)]);
Uy = zeros([1,length(a)]);
fx_force = zeros([length(a), N]);
fy_force = zeros([length(a), N]);

for J=1:length(a)
    a(J)
   
%%% prescribed boundary condition (tangential velocity) on the surface of squirmer:              

    V_surface = B1 * sin(angles) + B2 * sin(2 * angles);

       
%%% we use Levine and McKintosh response function
%%% v_\alpha = \alpha_{\alpha\beta} f_\beta
%%% in the code the response function is non-dimensionalized by multiplying
%%% it by 4\pi\eta h, where \eta is (3D) viscosity of the membrane and h is
%%% the thickness of the membrane: \tilde\alpha_\alpha\beta = alpha_\alpha\beta * 4\pi \eta h 
%%% Thus, in the code we actually use v_\alpha = \tilde\alpha_{\alpha\beta}
%%% * \tilde{f}_\beta, with \tilde{f}_\beta= \equiv f_\beta/ (4\pi\eta h),
%%% and \tilde{f}_\beta has units of m/s.


%%% Working in the lab frame we compute the forces on the squirmer as well
%%% as the translational velocity Ux and Uy 
       [fx_force(J,:), fy_force(J,:), Ux(J), Uy(J)] = ...
               solve_U_continuous(angles, a(J), Del, V_surface);
       
       FxNet = sum(fx_force(J,:)) %%% x-component of net force on squirmer
       FyNet = sum(fy_force(J,:)) %%% y-component of net force on squirmer

   Ux
   Uy
     
end

%%% Plot swimming speed rescaled by B1/2 vs a\ell_s
figure(1)
semilogx(a, Ux/(B1/2), 'o', 'LineWidth', 3)
hold on


m = 1%length(a) - 10; %%% m-th element of a-array
span_r = linspace(0, 3*a(m), 20);
  
 vx = zeros([1,length(span_r)]);
 vy = zeros([1,length(span_r)]);
 x=vx;
 y=vx;
 phi=0;  %%% choose direction along the direction of net motion
for i=1:length(span_r) 
         x(i) = span_r(i) * cos(phi);
         y(i) = span_r(i) * sin(phi);
         vx(i) = VX_FIELD(fx_force(m,:), fy_force(m,:), x(i), y(i), a(m));
         vy(i) = VY_FIELD(fx_force(m,:), fy_force(m,:), x(i), y(i), a(m));     
end

figure(2)
plot(span_r/a(m), -vx/(B1/2), 'ro', 'LineWidth', 3)
hold on
plot(span_r/a(m), vy, 'bo', 'LineWidth', 3)
hold off

% %%% Plot vectors of v-field at circumference of the circle
%  m = 1 %length(a) - 10; %%% m-th element of a-array
%  span_angles = linspace(0, 2*pi, 20);
%  x = a(m) * cos(span_angles);
%  y = a(m) * sin(span_angles);
%  
%  vx = zeros([1, length(x)]);
%  vy = zeros([1, length(y)]);
%  
%  for i=1:length(span_angles)
%      vx(i) = VX_FIELD(fx_force(m,:), fy_force(m, :), x(i), y(i), a(m));
%      vy(i) = VY_FIELD(fx_force(m,:), fy_force(m,:), x(i), y(i), a(m));     
%  end
%  
%  figure(4)
%  quiver(x, y, vx, vy, 'b')
%  hold off
%%%________________________________________________________________ 
%  span_r = linspace(1.5*a(1), 4*a(1), 10);
%  span_angles = linspace(0, 2*pi - 2*pi/10, 10);
% 
%  
%  vx = zeros([length(span_r), length(span_angles)]);
%  vy = zeros([length(span_r), length(span_angles)]);
%  x=vx;
%  y=vx;
%  
% for i=1:length(span_r) 
%     for j=1:length(span_angles)
%          x(i,j) = span_r(i) * cos(span_angles(j));
%          y(i,j) = span_r(i) * sin(span_angles(j));
%          vx(i,j) = VX_FIELD(fx, fy, x(i,j), y(i,j), a(1));
%          vy(i,j) = VY_FIELD(fx, fy, x(i,j), y(i,j), a(1));     
%     end
% end
% 
% 
%  figure(3)
%  quiver(x, y, vx, vy, 'b')
%  hold off
%%%___________________________________________________________________
%%% Case 3: plot  vector field  
 
m = 1 %length(a) - 10; %%% m-th element of a-array
 
x = [-5 * a(m) :  0.25 * a(m) : 5 * a(m)]';  %%% make a column
y =  -5 * a(m) :  0.25 * a(m) : 5 * a(m);    %%% make a row

X = repmat(x, [1,length(y)]);   %% form a matrix 
Y = repmat(y, [length(x), 1]);  %% form a matrix

fx_force(m,:)
fy_force(m,:)

VX = VX_FIELD(fx_force(m,:), fy_force(m,:), x, y, a(m));
VY = VY_FIELD(fx_force(m,:), fy_force(m,:), x, y, a(m));

figure(3)
 rectangle('Position',[-a(m), -a(m), 2*a(m), 2*a(m)],...
           'Curvature',[1,1],...
           'LineWidth',2,'LineStyle', '--', 'EdgeColor', 'r')
 daspect([1,1,1])
 hold on
 quiver(X, Y, VX, VY, 'b')
 hold off



end