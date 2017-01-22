% Plots the trajectory for a boundary element squirmer in an enclosed domain 

tic


%These are scaled to become the default simulation time T and timestep
%dt_o, although both will change depending on how dt is rescaled during the
%simulation.
Time = 70;
Increment = 1;



% a_array = [1 0.01 0.0001]; %[1e-3 1e-2 1 10];
% R_array = [10]; %[5 10 25]; This is the ratio R/a
% B2_array = [-1 0 1]; %[-5s -1 0 1 5];
% B1_array = [1]; %[1 0];



%When h ~=0 an if statement replaces the circular enclosure with a vertical
%wall

a_array = [0.1];
R_array = [10];
B2_array = [-1];
B1_array = [0];

h = 15*a; % h is the height of a vertical wall
%that is placed symmetrically about the origin at a distance R to the
%right of the origin

FinishedLoopCount = 0;

for i = 1:length(a_array) %Cycle through beast radii
    for j = 1:length(R_array) %Cycle through enclosure radii
        for k = 1:length(B2_array) %Cycle B2
            for m = 1:length(B1_array) %Cycle B1
                close all
                

                
                a = a_array(i);  %%% radius of the disk nondimensionalized by the Saffman length
                R = R_array(j)*a;  %%% Radius of enclosure nondimensionalized by the Saffman length
                B2 = B2_array(k);
                B1 = B1_array(m);
                
                s = 0.1 * a;          %%% radial spacing between neighboring blobs
                epsilon = s/8;        %%% radius of blobs
                
                %%%Needs to be implemented into solve_U_enclosure
                Scale = 10/a;         %%% This scales some terms in the matrix equation. Possibly helps with invertibility.
            

                %This rescaling makes the trajectories similar in distance
                %when compared to the enclosure size. (This way I can cycle
                %through many beast radii and leave [Time] fixed.)
                T = Time*a;         
                dt = Increment*a;  

                %Initial Conditions
                r_o = .3*R;         %%% Radial coordinate of beast cm from center of enclosure
                phi_o = 0;       %%% Angle coordinate of beast cm from center of enclosure
                theta_o = pi/3;   %%% Beast intial orientation (head direction in enclosure frame)

                x_o = r_o*cos(phi_o); %%% Beast CM initial x position as seen in enclosure frame.
                y_o = r_o*sin(phi_o); %%% Beast CM Initial y position in enclosure frame.



                %Coordinates of beast blobs in beast frame.
                %This places a beast in the center of the enclosure.
                [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

                Nbeast = sum(BlobsPerLayer); %%% Number of blobs in the beast
                NRim = BlobsPerLayer(end);   %%% Number of blobs in the outermost beast layer

                %Rotate coordinates so beast faces according to theta_o
                [xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, theta_o);

                %Prescribe tangential velocities
                %And Rotate velocities according to theta_o into enc frame.
                [VxRim, VyRim] = UpdatedPrescribeWave(NRim, B1, B2, theta_o);

                %Enclosure Blob Coordinates.
                %Enclosure is a ring centered about the origin.
                
                % Primary Enclosure Discretization
                d = s; %Enclosure blob spacing same as beast blob spacing
                [x_Enc, y_Enc] = DiscretizeEnclosure(R,d);
                %
                
                %Replace enclosure coordinates with wall if selected
                if h ~= 0
                    [x_Enc, y_Enc] = DiscretizeWall(R,s,h);
                end
                
                NEnc = length(y_Enc);
                
                %Alternate enclosure discretization. Allows adjustments of
                %the number of enclosure blobs. (Thought this might help
                %with collisions.)
%                 Nencl = 500;
%                 x_Enc = R * cos(linspace(0, 2*pi, Nencl));
%                 y_Enc = R * sin(linspace(0, 2*pi, Nencl));

                %Translate beast geometric center to the chosen initial position in enclosure frame.
                xcoord = xcoord + x_o;
                ycoord = ycoord + y_o;

                %Grab the head coordinate for plotting
                x_head = xcoord(end - NRim + 1);
                y_head = ycoord(end - NRim + 1);

                
                %Begin Simulation
                %Right now it returns way too much
                [Ux_history, Uy_history, W_history, theta_history, x_cm_history, y_cm_history, separation_history, ...
                    dt_history, time_history, x_history, y_history, Matrix_history,x_cm_history_Recip, y_cm_history_Recip, ...
                    separation_history_Recip,speed_history_Recip, W_history_Recip, x_history_Recip, y_history_Recip]...
                    = TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, R, a, Scale, B1, B2);

                %For checking the speed relative to distance from edge.
                speed_history = (Uy_history.^2 + Ux_history.^2).^(1/2);

                
                %Display calculation percentage
                FinishedLoopCount = FinishedLoopCount + 1;
                LoopPercentage = FinishedLoopCount/(length(a_array)*length(R_array)*length(B1_array)*length(B2_array))*100

                toc

                %% Create some strings for plot detail
                str_T = ['Total Time = ',num2str(T)];
                str_dt = ['Base dt = ',num2str(dt)];
                str_Time = 'Nondimen by B_1/(2*l_s)';
                str_a = ['a = ',num2str(a)];
                str_s = ['s = ',num2str(s)];
                str_eps = ['epsilon = ',num2str(epsilon)];
                str_R = ['R = ',num2str(R)];
                str_r_o = ['r_o = ',num2str(r_o)];
                str_phi_o = ['phi_o = ',num2str(phi_o)];
                str_theta_o = ['theta_o = ',num2str(theta_o)];
                str_B1 = ['B1 = ',num2str(B1)];
                str_B2 = ['B2 = ',num2str(B2)];

                descr = {'Parameters:';
                    str_T;
                    str_dt;
                    str_Time;
                    str_a;
                    str_s;
                    str_eps;
                    str_R;
                    str_r_o;
                    str_phi_o;
                    str_theta_o;
                    str_B1;
                    str_B2};
                
                descr = {'Parameters:';
                    str_B1;
                    str_B2};

                %% Plot the cm trajectory
                fig = figure(1);
                ax1 = axes('Position',[0 0 1 1],'Visible','off');
                ax2 = axes('Position',[.3 .1 .6 .8]);
                axes(ax1);
                text(.1,0.6,descr);
                %back to data
                axes(ax2);

                scatter(x_cm_history(2:end-2), y_cm_history(2:end-2), '.', 'k');
                hold on
%                 scatter(x_cm_history_Recip(2:end-2), y_cm_history_Recip(2:end-2), '.', 'm');
%                 legend('Boundary Element', 'Lorentz Reciprocal','Location','southwest');
                scatter(x_cm_history(1), y_cm_history(1), '.', 'g'); %Begin at green
                scatter(x_cm_history(end-1), y_cm_history(end-1), '.', 'y'); %End at yellow
                
                %Body. Inner then rim. At two points in time (1 and end)
                scatter(x_history(end,1:Nbeast-NRim),y_history(end,1:Nbeast-NRim),'.','b')
                scatter(x_history(end,Nbeast-NRim+1:end),y_history(end,Nbeast-NRim+1:end),'.','r')
%                
                scatter(x_history(1,1:Nbeast-NRim),y_history(1,1:Nbeast-NRim),'.','b')
                scatter(x_history(1,Nbeast-NRim+1:end),y_history(1,Nbeast-NRim+1:end),'.','r')
                
                scatter(x_history(floor(end/2),1:Nbeast-NRim),y_history(floor(end/2),1:Nbeast-NRim),'.','b')
                scatter(x_history(floor(end/2),Nbeast-NRim+1:end),y_history(floor(end/2),Nbeast-NRim+1:end),'.','r')
                
                %Head
                scatter(x_history(end,Nbeast - NRim + 1),y_history(end,Nbeast - NRim + 1),400,'.','m')
                scatter(x_history(1,Nbeast - NRim + 1),y_history(1,Nbeast - NRim + 1),400,'.','m')
                scatter(x_history(floor(end/2),Nbeast - NRim + 1),y_history(floor(end/2),Nbeast - NRim + 1),400,'.','m')
                
                %Lorentz data
                
%                  scatter(x_cm_history_Recip(1), y_cm_history_Recip(1), '.', 'g'); %Begin at green
%                  scatter(x_cm_history_Recip(end-1), y_cm_history_Recip(end-1), '.', 'r'); %End at red
                
                %Body
%                 scatter(x_history_Recip(3,:),y_history_Recip(3,:),'.')
%                 scatter(x_history_Recip(1,:),y_history_Recip(1,:),'.')
                %Head
%                 scatter(x_history_Recip(3,Nbeast - NRim + 1),y_history_Recip(3,Nbeast - NRim + 1),'.')
%                 scatter(x_history_Recip(1,Nbeast - NRim + 1),y_history_Recip(1,Nbeast - NRim + 1),'.')
                

                
                
                
                daspect([1,1,1])
                scatter(x_Enc, y_Enc, '.');
                axis off
                hold off
                title1 = ['trajectory' num2str(i) num2str(j) num2str(k) num2str(m)];
                saveas(gcf,title1, 'png')
                saveas(gcf,title1, 'eps')
                %% Plotting velocity vs separation
%                 fig = figure(2);
% %                 ax1 = axes('Position',[0 0 1 1],'Visible','off');
% %                 ax2 = axes('Position',[.3 .1 .6 .8]);
% %                 axes(ax1);
% %                 text(.025,0.6,descr);
% %                 %back to plotting data
% %                 axes(ax2);
%                 hold on
%                 scatter(separation_history(2:end-1), speed_history(2:end-1));
%                 scatter(separation_history_Recip(2:end-1), speed_history_Recip(2:end-1));
%                 legend('Boundary Element', 'Lorentz Reciprocal', 'Location','northeast');
%                 xlabel('Nondimensional Separation Distance [a/l_s]')
%                 ylabel('Nondimensional Swimming Speed [U/(B_1/2)]')
%                 hold off
%                 title2 = ['speed' num2str(i) num2str(j) num2str(k) num2str(m)];
%                 saveas(gcf,title2, 'png')
%                 saveas(gcf,title2, 'eps')
                %% Plotting angular velocity vs separation
%                 fig = figure(3);
% %                 ax1 = axes('Position',[0 0 1 1],'Visible','off');
% %                 ax2 = axes('Position',[.3 .1 .6 .8]);
% %                 axes(ax1);
% %                 text(.025,0.6,descr);
% %                 %back to plotting data
% %                 axes(ax2);
%                 hold on
%                 scatter(separation_history(2:end-1), W_history(2:end-1));
%                 scatter(separation_history_Recip(2:end-1), W_history_Recip(2:end-1));
%                 legend('Boundary Element', 'Lorentz Reciprocal', 'Location','northeast');
%                 xlabel('Nondimensional Separation Distance (a/l_s)')
%                 ylabel('Nondimensional Angular Velocity [\Omega /(B_1/2)]')
%                 hold off
%                 title3 = ['angular' num2str(i) num2str(j) num2str(k) num2str(m)];
%                 saveas(gcf,title3, 'png')
%                 saveas(gcf,title3, 'eps')
                %% Plotting separation vs simulation time
%                 fig = figure(4);
% %                 ax1 = axes('Position',[0 0 1 1],'Visible','off');
% %                 ax2 = axes('Position',[.3 .1 .6 .8]);
% %                 axes(ax1);
% %                 text(.025,0.6,descr);
% %                 %back to plotting data
% %                 axes(ax2);
%                 scatter(time_history(2:end-1), separation_history(2:end-1));
%                 %title('Distance from Enclosure Wall vs Time','FontSize',16,'FontWeight','bold')
%                 xlabel('Nondimensional Time [(B_1/(2 l_s) t]')
%                 ylabel('Nondimensional Separation [a/l_s]')
%                 title4 = ['separation' num2str(i) num2str(j) num2str(k) num2str(m)];
%                 saveas(gcf,title4, 'png')
%                 saveas(gcf,title4, 'eps')
            end
        end
    end
end