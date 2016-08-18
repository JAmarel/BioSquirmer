% Second inactive attempt. Pick a U and solve for F by summing the little f's.
% Active and Inactive disks must have velocity in the same direction (x)

% 1 
tic
n = 20; %Radius sampling
Radii = logspace(-1.5,3,n); % logspace(a,b,n) generates n points between decades 10^a and 10^b.

m = 1; %Spacing sampling


ActiveUx = zeros([n, m]);
ActiveUy = zeros([n,m]);
ActiveU = zeros([n,m]);

RecipUx = zeros([n, m]);
RecipUy = zeros([n,m]);
RecipU = zeros([n,m]);

for i=1:n
    a = Radii(i);
        Spacings = linspace(.1, .1, m);
    for j=1:m
        s = Spacings(j)*a;         %%% spacing between neighboring blobs
        epsilon = s/8;      %%% radius of the blob
        theta_o = 0.01;       %%% Beast intial orientation (head direction)

        %Velocity of inactive disk
        U_hat = 1;
        
        Ux_hat = U_hat*cos(theta_o); 
        Uy_hat = U_hat*sin(theta_o);

        [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

        Nblobs = sum(BlobsPerLayer); %%% total number of blobs 

        NR = length(BlobsPerLayer); %%% Number of radial layers
        NRim = BlobsPerLayer(end);  %%% number of blobs in the outermost layer

        [VxRim, VyRim, B1] = PrescribeWave(NRim);
        
        [VxRim, VyRim] = Rotate_Vector(VxRim, VyRim, theta_o);

        %Nondimensionalize
        VxRim = VxRim/(B1/2);
        VyRim = VyRim/(B1/2);

        %% Inactive Recip Thm
        
        % Hat is conjugate scenario
        % Should I use the equations below, or FxNet_hat?
        Fx_hat = Ux_hat/HPW_mobility(a);
        Fy_hat = Uy_hat/HPW_mobility(a);
        
        % First, solve for Ux as predicted by recip thm
        [fx_hat_1, fy_hat_1] = solve_U_disk_inactive(xcoord, ycoord, epsilon, NRim, Ux_hat, 0);

        FxRim_hat_1 = fx_hat_1(end-NRim+1:end);
        FyRim_hat_1 = fy_hat_1(end-NRim+1:end);

        FxNet_hat_1 = sum(fx_hat_1); %%% x-component of net force on squirmer
        FyNet_hat_1 = sum(fy_hat_1); %%% y-component of net force on squirmer
        
        Ux_Recip = (-1/(Fx_hat))*(sum(VxRim.*FxRim_hat_1) + sum(VyRim.*FyRim_hat_1));
        
        %Ux_Recip = (-1/(FxNet_hat_1))*(sum(VxRim.*FxRim_hat_1) + sum(VyRim.*FyRim_hat_1)); % This line is probably incorrect. (Gives perfect agreement though.)
        
        % Second, solve for Uy from recip thm
        [fx_hat_2, fy_hat_2] = solve_U_disk_inactive(xcoord, ycoord, epsilon, NRim, 0, Uy_hat);
        
        FxRim_hat_2 = fx_hat_2(end-NRim+1:end);
        FyRim_hat_2 = fy_hat_2(end-NRim+1:end);

        FxNet_hat_2 = sum(fx_hat_2); %%% x-component of net force on squirmer
        FyNet_hat_2 = sum(fy_hat_2); %%% y-component of net force on squirmer
        
        Uy_Recip = (-1/(Fy_hat))*(sum(VxRim.*FxRim_hat_2) + sum(VyRim.*FyRim_hat_2));
        
        %Uy_Recip = (-1/(FyNet_hat_2))*(sum(VxRim.*FxRim_hat_2) + sum(VyRim.*FyRim_hat_2)); %Possibly incorrect
        
        % Third, solve for W from recip thm.
        %% Active
        % Solve for Ux as predicted by active swimmer code
        
        [fx, fy, Ux, Uy, W, Matrix] = solve_U_disk_rot(xcoord, ycoord, epsilon, VxRim, VyRim, NRim);

        FxRim = fx(end-NRim+1:end);
        FyRim = fy(end-NRim+1:end);

        FxNet = sum(fx); %%% x-component of net force on squirmer
        FyNet = sum(fy); %%% y-component of net force on squirmer

        TorqueNet2 = dot(xcoord,fy.') - dot(ycoord,fx.');

        
        %% Fill in data
        ActiveUx(i,j) = Ux;
        ActiveUy(i,j) = Uy;
        ActiveU(i,j) = sqrt(Ux^2 + Uy^2);
        
        % If statement is necessary to avoid dividing by zero
        if Uy_hat == 0
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
    end
end
toc

        %% Plotting
for j=1:m
    fig = figure(j);
    str_spacing = ['Spacing (s/a) = ', num2str(Spacings(j))];
    descr = {'Parameters:';
    str_spacing};
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    ax2 = axes('Position',[.3 .1 .6 .8]);
    axes(ax1);
    text(.025,0.6,descr);
    %back to data
    axes(ax2);
    semilogx(Radii, ActiveU(:,j), 'k-')
    hold on
    semilogx(Radii, RecipU(:,j), 'bo')
    title('Swimming Speed vs. Radius','FontSize',16,'FontWeight','bold')
    xlabel('Log Scale Nondimensional Radius (a/l_s)')
    ylabel('Nondimensional Swimming Speed U/(B1/2)')
    hold off
    %saveas(gcf,[num2str(j),'.png']);
end