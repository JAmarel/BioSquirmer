%Looking for swimming speed dependence on number of rim blobs/packing density

%DiscretizeDisk_Pack needs slight changes for this script to function
%again.

n = 20;

a = 10;
s = 0.08 * a;          %%% spacing between neighboring blobs
epsilon = s/8;       %%% radius of the blob
B1 = 1;
B2 = 0;

CircumSpacing = linspace(s/(3),1.2*s,n); %slight overlap if d < e/2

Efficiencies = zeros([1, n]);
Blobs = zeros([1,n]);
Speeds = zeros([1, n]);

for i=1:n
        %%% spacing between neighboring blobs
    d = CircumSpacing(i);          %%% blob spacing
      
    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk_Pack(a,s,d);
    
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
    Blobs(i) = Nblobs;
    CircumSpacing(i) = d;
    Speeds(i) = speed;
    
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
figure(1)
plot(Blobs, Speeds, 'o')
title('Swimming Speed vs. Total Number of Blobs','FontSize',16,'FontWeight','bold')
xlabel('Number of Blobs in Beast')
ylabel('Nondimensional Speed')
%saveas(gcf,'Efficiency vs Nrim.png')

% figure(2)
% plot(CircumSpacing, Speeds, 'o')
% title('Swimming Speed vs. Blob Spacing','FontSize',16,'FontWeight','bold')
% xlabel('Circumferential Blob Spacing Nondimensional [a/l_s]')
% ylabel('Nondimensional Speed')
% %saveas(gcf,'Efficiency vs d.png')

% figure(3)
% loglog(Radii, Efficiencies, 'o')
% title('Log-Log Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
% xlabel('Log Scale Nondimensional Radius (a/l_s)')
% ylabel('Log Scale Efficiency (P_d_r_a_g/P_s_w_i_m)')
% saveas(gcf,'Efficiency 3.png')
% ylabel('Log Scale Efficiency (P_d_r_a_g/P_s_w_i_m)')


