%Looking for swimming efficiency dependence on number of rim blobs/packing density

n = 40;

a = 10;
s = 0.08 * a;          %%% spacing between neighboring blobs
epsilon = s/8;       %%% radius of the blob

CircumSpacing = linspace(s/(3),1.2*s,n); %slight overlap if d < e/2

Efficiencies = zeros([1, n]);
RimBlobs = zeros([1,n]);

for i=1:n
        %%% spacing between neighboring blobs
    d = CircumSpacing(i);          %%% blob spacing
      
    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk_Pack(a,s,d);
    
    Nblobs = sum(BlobsPerLayer); %%% total number of blobs 
    NR = length(BlobsPerLayer); %%% Number of radial layers
    NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer
    
    [VxRim, VyRim, B1] = PrescribeWave(NRim);
    
    
    [fx, fy, Ux, Uy] = solve_U_disk(xcoord, ycoord, epsilon, VxRim, VyRim, NRim); %Currently not calling the rotating solver.
    
    fx = fx/(B1/2); %Nondimensionalizing.
    fy = fy/(B1/2);
    Ux = Ux/(B1/2);
    Uy = Uy/(B1/2);
    VxRim = VxRim/(B1/2);
    VyRim = VyRim/(B1/2);

    FxRim = fx(end-NRim+1:end);
    FyRim = fy(end-NRim+1:end);

    FxNet = sum(fx); %%% x-component of net force on squirmer
    FyNet = sum(fy); %%% y-component of net force on squirmer

    speed = sqrt(Ux^2 + Uy^2);

    efficiency = CalcEfficiency(FxRim, FyRim, VxRim, VyRim, a, speed)*(B1/2)^2;
    Efficiencies(i) = efficiency;
    RimBlobs(i) = NRim;
    CircumSpacing(i) = d;
    
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
plot(RimBlobs, Efficiencies, 'o')
title('Swimming Efficiency vs. Number of Rim Blobs','FontSize',16,'FontWeight','bold')
xlabel('Number of Rim Blobs')
ylabel('Efficiency (P_d_r_a_g/P_s_w_i_m)')
%saveas(gcf,'Efficiency vs Nrim.png')

figure(2)
plot(CircumSpacing, Efficiencies, 'o')
title('Swimming Efficiency vs. Blob Spacing','FontSize',16,'FontWeight','bold')
xlabel('Circumferential Blob Spacing Nondimensional [a/l_s]')
ylabel('Efficiency (P_d_r_a_g/P_s_w_i_m)')
%saveas(gcf,'Efficiency vs d.png')

% figure(3)
% loglog(Radii, Efficiencies, 'o')
% title('Log-Log Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
% xlabel('Log Scale Nondimensional Radius (a/l_s)')
% ylabel('Log Scale Efficiency (P_d_r_a_g/P_s_w_i_m)')
% saveas(gcf,'Efficiency 3.png')
% ylabel('Log Scale Efficiency (P_d_r_a_g/P_s_w_i_m)')

