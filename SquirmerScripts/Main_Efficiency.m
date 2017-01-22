%%Calculates/plots the swimming efficiency for various beast radii
n = 30; %Number of radii to sample
Radii = logspace(-8,8,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

Efficiencies = zeros([1, n]);
ShellNumber = zeros([1, n]);

B1 = 1;
B2 = 0;

for i=1:n
    a = Radii(i);
    s = 0.08 * a;          %%% spacing between neighboring blobs
    epsilon = s/8;       %%% radius of the blob
    
    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);
    
    Nblobs = sum(BlobsPerLayer); %%% total number of blobs 
    NR = length(BlobsPerLayer); %%% Number of radial layers
    NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer
    
    [VxRim, VyRim] = PrescribeWave(NRim, B1, B2);
    
    %Nondimensionalize
    VxRim = VxRim/(B1/2);
    VyRim = VyRim/(B1/2);
    
    [fx, fy, Ux, Uy] = solve_U_disk(xcoord, ycoord, epsilon, VxRim, VyRim, NRim); %Currently not calling the rotating solver.
    

    FxRim = fx(end-NRim+1:end);
    FyRim = fy(end-NRim+1:end);

    FxNet = sum(fx); %%% x-component of net force on squirmer
    FyNet = sum(fy); %%% y-component of net force on squirmer

    speed = sqrt(Ux^2 + Uy^2);

    efficiency = CalcEfficiency(FxRim, FyRim, VxRim, VyRim, a, speed)*(B1/2)^2;
    Efficiencies(i) = efficiency;
    ShellNumber(i) = NR;
    
%%Plot the position blobs by xcoord and ycoord
% figure(i+3)
% plot(xcoord, ycoord, 'o')
% str1 = strcat('a = ',num2str(a));
% str2 = strcat('NR = ',num2str(NR));
% text(1.1*max(xcoord),1.1*max(ycoord),str1,'FontSize',20);
% text(-1.1*max(xcoord),1.1*max(ycoord),str2,'FontSize',20);
% daspect([1,1,1])
% hold on
% plot(xcoord(Nblobs - NRim + 1:end), ycoord(Nblobs - NRim + 1:end), 'ro', 'LineWidth', 3)
% axis off
% hold off
end

%%Plot the efficiency vs nondimensional radius
% figure(1)
% plot(Radii, Efficiencies, 'o')
% title('Swimming Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
% xlabel('Nondimensional Radius (a/l_s)')
% ylabel('Efficiency (P_d_r_a_g/P_s_w_i_m)')
% saveas(gcf,'Efficiency 1.png')

figure(2)
semilogx(Radii, Efficiencies, 'o')
title('Swimming Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
xlabel('Nondimensional Beast Radius [a/l_s]')
ylabel('Efficiency (P_d_r_a_g/P_s_w_i_m)')
saveas(gcf,'EfficiencyvsRadius.png')
saveas(gcf,'EfficiencyvsRadius.eps')

% figure(3)
% loglog(Radii, Efficiencies, 'o')
% title('Log-Log Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
% xlabel('Log Scale Nondimensional Radius (a/l_s)')
% ylabel('Log Scale Efficiency (P_d_r_a_g/P_s_w_i_m)')
% saveas(gcf,'Efficiency 3.png')

% figure(4)
% plot(ShellNumber, Efficiencies, 'o')
% x = [0.5 0.5];
% y = [0.3 0.6];
% annotation('textarrow',x,y,'String','Increasing Beast Radius ','FontSize',7)
% xlim([10 15])
% title('Swimming Efficiency vs Number of Blob Shells','FontSize',16,'FontWeight','bold')
% xlabel('Number of Blob Shells')
% ylabel('Efficiency (P_d_r_a_g/P_s_w_i_m)')
% saveas(gcf,'ShellRound.png')


