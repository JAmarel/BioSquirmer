n = 20;
Radii = logspace(-2,2,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

Efficiencies = zeros([1, n]);

for i=1:n
    a = Radii(i);
    s= 0.1 * a;          %%% spacing between neighboring blobs
    epsilon = s/8;       %%% radius of the blob
    
    [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);
    Nblobs = sum(BlobsPerLayer); %%% total number of blobs 
    
    NR = length(BlobsPerLayer); %%% Number of radial layers
    NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer
    
    [VxRim, VyRim, B1] = PrescribeWave(NRim);
    [fx, fy, Ux, Uy] = solve_U_disk(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);

    FxRim = fx(end-NRim+1:end);
    FyRim = fy(end-NRim+1:end);

    FxNet = sum(fx); %%% x-component of net force on squirmer
    FyNet = sum(fy); %%% y-component of net force on squirmer

    speed = Ux/(B1/2);   %%%% swimming velocity non-dimensionalized by B1/2

    efficiency = CalcEfficiency(FxRim, FyRim, VxRim, VyRim, a, speed)*(B1/2)^2;
    Efficiencies(i) = efficiency;
end

%%Plot the efficiency vs nondimensional radius
figure(1)
plot(Radii, Efficiencies, 'o')
title('Swimming Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
xlabel('Nondimensional Radius (a/l_s)')
ylabel('Efficiency (P_d_r_a_g/P_s_w_i_m)')
saveas(gcf,'Efficiency 1.png')

figure(2)
semilogx(Radii, Efficiencies, 'o')
title('Swimming Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
xlabel('Log Scale Nondimensional Radius (a/l_s)')
ylabel('Efficiency (P_d_r_a_g/P_s_w_i_m)')
saveas(gcf,'Efficiency 2.png')

figure(3)
loglog(Radii, Efficiencies, 'o')
title('Log-Log Efficiency vs. Radius','FontSize',16,'FontWeight','bold')
xlabel('Log Scale Nondimensional Radius (a/l_s)')
ylabel('Log Scale Efficiency (P_d_r_a_g/P_s_w_i_m)')
saveas(gcf,'Efficiency 3.png')
ylabel('Log Scale Efficiency (P_d_r_a_g/P_s_w_i_m)')
