% Computes and compares the swimming speed vs radius between the boundary element approach and the numerical reciprocal theorem for an unbounded domain.
% loops through variables n and m, to cycle beast radii and blob spacing,
% could be modified to cycle B1/B2

% Active and Inactive disks must have velocity in the same direction?

% 1 
tic
n = 20; %Number of radii to sample
Radii = logspace(-1.5,3,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

B1 = 1;
B2 = 0;

m = 1; %Number of spacings to sample

HPWMobility = zeros([n,m]);
NumericalMobility = zeros([n,m]);

ActiveUx = zeros([n, m]);
ActiveUy = zeros([n,m]);
ActiveU = zeros([n,m]);

RecipUx = zeros([n, m]);
RecipUy = zeros([n,m]);
RecipU = zeros([n,m]);

for i=1:n %Loop through non-dim radii
    a = Radii(i);
        Spacings = linspace(.1, .1, m);
    for j=1:m %Loop through blob spacing
        s = Spacings(j)*a;         %%% spacing between neighboring blobs
        epsilon = s/8;      %%% radius of a blob
        theta_o = 0;       %%% Beast intial orientation (head direction)

        %Velocity of inactive disk
        %Choose simplicity
        U_hat = 1;
        W_hat = 1;
        
        
        Ux_hat = U_hat*cos(theta_o); 
        Uy_hat = U_hat*sin(theta_o);

        [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s); %Initiate a beast at the origin
        
        %Rotate coordinates according to theta_o
        [xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, theta_o);

        Nblobs = sum(BlobsPerLayer); %%% total number of blobs 

        NR = length(BlobsPerLayer); %%% Number of radial layers
        NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

        %Prescribe wave, nondimensionalize, and rotation according to theta_o.
        [VxRim, VyRim] = UpdatedPrescribeWave(NRim, B1, B2, theta_o);
       

        %% Inactive Recip Thm
        
        % Hat is conjugate scenario
        % Should I use the equations below, or FxNet_hat_1?
        % This choice is labeled alternate calculation
        
        % The value for Fx_hat is the resulting net force on the body necessary for our
        % choice of Ux_hat and Uy_hat according to the predictions of our
        % HPW mobility calculation.
        Fx_hat = Ux_hat/HPW_mobility(a);
        Fy_hat = Uy_hat/HPW_mobility(a);
        
        %The alternate choice is a similar object, the resulting net force
        %on the body necessary for our choice of Ux_hat and Uy_hat, but
        %this time according to the results of the Lorentz Reciprocal
        %theorem solver. I chose this option for now, because as we have to
        %use the rim forces from the lorentz solver, it seems best for
        %consistency to also use the total force of the same origin.%The
        %values are both reasonable and similar.
        
        
        
        % Recip theorum allows us to solve for Ux, Uy, and W separately.
        
        % First, solve for Ux as predicted by recip thm
        [fx_hat_1, fy_hat_1] = solve_U_disk_inactive(xcoord, ycoord, epsilon, NRim, Ux_hat, 0, 0);

        FxRim_hat_1 = fx_hat_1(end-NRim+1:end);
        FyRim_hat_1 = fy_hat_1(end-NRim+1:end);

        FxNet_hat_1 = sum(fx_hat_1); %%% x-component of net rim force on squirmer
        FyNet_hat_1 = sum(fy_hat_1); %%% y-component of net rim force on squirmer
        
        %Alternate
        %Ux_Recip = (-1/(Fx_hat))*(sum(VxRim.*FxRim_hat_1) + sum(VyRim.*FyRim_hat_1));
        
        Ux_Recip = (-1/(FxNet_hat_1))*(sum(VxRim.*FxRim_hat_1) + sum(VyRim.*FyRim_hat_1));
        
        % Second, solve for Uy from recip thm
        [fx_hat_2, fy_hat_2] = solve_U_disk_inactive(xcoord, ycoord, epsilon, NRim, 0, Uy_hat, 0);
        
        FxRim_hat_2 = fx_hat_2(end-NRim+1:end);
        FyRim_hat_2 = fy_hat_2(end-NRim+1:end);

        FxNet_hat_2 = sum(fx_hat_2); %%% x-component of net rim force on squirmer
        FyNet_hat_2 = sum(fy_hat_2); %%% y-component of net rim force on squirmer
        
        %Alternate
        %Uy_Recip = (-1/(Fy_hat))*(sum(VxRim.*FxRim_hat_2) + sum(VyRim.*FyRim_hat_2));
        
        Uy_Recip = (-1/(FyNet_hat_2))*(sum(VxRim.*FxRim_hat_2) + sum(VyRim.*FyRim_hat_2));
        
        % Third, solve for W from recip thm.

  
        [fx_hat_3, fy_hat_3] = solve_U_disk_inactive(xcoord, ycoord, epsilon, NRim, 0, 0, W_hat);
        
        FxRim_hat_3 = fx_hat_3(end-NRim+1:end);
        FyRim_hat_3 = fy_hat_3(end-NRim+1:end);
        
        Torque_hat = dot(xcoord,fy_hat_3.') - dot(ycoord,fx_hat_3.');
        
        W_Recip = (-1/(Torque_hat))*(sum(VxRim.*FxRim_hat_3) + sum(VyRim.*FyRim_hat_3));
        %% Active
        % Solve for Ux as predicted by active swimmer code
        
        [fx, fy, Ux, Uy, W, Matrix] = solve_U_disk_rot(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);

        FxRim = fx(end-NRim+1:end);
        FyRim = fy(end-NRim+1:end);

        FxNet = sum(fx); %%% x-component of net force on squirmer
        FyNet = sum(fy); %%% y-component of net force on squirmer

        TorqueNet = dot(xcoord,fy.') - dot(ycoord,fx.');

        
        %% Fill in data
        ActiveUx(i,j) = Ux;
        ActiveUy(i,j) = Uy;
        ActiveU(i,j) = sqrt(Ux^2 + Uy^2);
        
        % If statement is necessary to avoid dividing by zero
        if Uy_hat == 0 %When Uy_hat=0, Fy_hat=0 and calculating Uy_Recip involves division by Fy_hat
            RecipUy(i,j) = 0;
        else
            RecipUy(i,j) = Uy_Recip;
        end
        
        if Ux_hat == 0
            RecipUx(i,j) = 0;
        else
            RecipUx(i,j) = Ux_Recip;
        end
        
        RecipU(i,j) = sqrt(RecipUy(i,j)^2 + RecipUx(i,j)^2);
        
        HPWMobility(i,j) = HPW_mobility(a);
        NumericalMobility(i,j) = U_hat/sqrt(FxNet_hat_1^2 + FyNet_hat_2^2);
    end
end
toc

        %% Plotting Speeds vs Beast Radius. Comparing
for j=1:m
    fig = figure(j);
%     str_spacing = 'Spacing (s/a):';
%     str_s = num2str(Spacings(j));
%     descr = {'Parameters:';
%     str_spacing;
%     str_s};
%     ax1 = axes('Position',[0 0 1 1],'Visible','off');
%     ax2 = axes('Position',[.3 .1 .6 .8]);
%     axes(ax1);
%     text(.025,0.6,descr);
%     %back to data
%     axes(ax2);
    semilogx(Radii, ActiveU(:,j), 'k-')
    hold on
    semilogx(Radii, RecipU(:,j), 'bo')
    title('Swimming Speed vs. Radius','FontSize',16,'FontWeight','bold')
    xlabel('Nondimensional Beast Radius [a/l_s]')
    ylabel('Nondimensional Swimming Speed [U/(B1/2)]')
    legend('Boundary Element Speed', 'Lorentz Speed');
    hold off
    %saveas(gcf,[num2str(j),'.png']);
    saveas(gcf,'LorentzSpeedvsRadius.png')
    saveas(gcf,'LorentzSpeedvsRadius.eps')
end

        %% Plotting Mobility vs Beast Radius. Comparing
% for j=1:m
%     fig = figure(j);
% %     str_spacing = ['Spacing (s/a) = ', num2str(Spacings(j))];
% %     descr = {'Parameters:';
% %     str_spacing};
% %     ax1 = axes('Position',[0 0 1 1],'Visible','off');
% %     ax2 = axes('Position',[.3 .1 .6 .8]);
% %     axes(ax1);
% %     text(.025,0.6,descr);
% %     back to data
% %     axes(ax2); %#ok<*LAXES>
%     scatter(Radii, HPWMobility(:,j), 'k', '.')
%     hold on
%     scatter(Radii, NumericalMobility(:,j), 'b')
%     set(gca,'xscale','log')
%     legend('HPW Mobility', 'Lorentz Mobility');
%     title('Approximated HPW Mobility vs. Numerical Mobility','FontSize',16,'FontWeight','bold')
%     xlabel('Nondimensional Beast Radius [a/l_s]')
%     ylabel('Reduced Mobility [4\pi \eta_m \mu]')
%     hold off
%     %saveas(gcf,num2str(j),'eps');
%     saveas(gcf,'MobilityvsRadius.png')
%     saveas(gcf,'MobilityvsRadius.eps')
% end