 tic

%These scale with a
Steps = 10;
Increment = .5;

a_array = [1e-3 1e-2 1 10];
R_array = [5 10 25];
B2_array = [-5 -1 0 1 5];
B1_array = [1 0];

FinishedLoopCount = 0;

for i = 1:length(a_array)
    for j = 1:length(R_array)
        for k = 1:length(B2_array)
            for m = 1:length(B1_array)
                close all
                

                
                a = a_array(i);  %%% radius of the disk nondimensionalized by the Saffman length
                R = R_array(j)*a;  %%% Radius of enclosure
                B2 = B2_array(k);
                B1 = B1_array(m);
                
                s = 0.1 * a;          %%% radial spacing between neighboring blobs
                epsilon = s/8;        %%% radius of blobs
                Scale = 10/a;
            


                T = Steps*a;
                dt = Increment*a;  

                %Initial Conditions
                r_o = .4*R;         %%% Radial coordinate of beast cm from center of enclosure
                phi_o = -.6*pi/2;       %%% Angle coordinate of beast cm from center of enclosure
                theta_o = pi/4;   %%% Beast intial orientation (head direction)

                x_o = r_o*cos(phi_o); %%% Beast CM initial x position as seen in enclosure frame.
                y_o = r_o*sin(phi_o); %%% Beast CM Initial y position



                %Coordinates of beast blobs in beast frame.
                [xcoord, ycoord, BlobsPerLayer] = DiscretizeDisk(a,s);

                Nbeast = sum(BlobsPerLayer); %%% Number of blobs in the beast
                NRim = BlobsPerLayer(end);   %%% Number of blobs in the outermost beast layer

                %Rotate coordinates according to theta_o
                [xcoord, ycoord] = Rotate_Vector(xcoord, ycoord, theta_o);

                %Prescribe wave in the beast frame.
                [VxRim, VyRim] = PrescribeWave(NRim, B1, B2);

                %Now Rotate velocities according to theta_o into lab frame.
                [VxRim, VyRim] = Rotate_Vector(VxRim, VyRim, theta_o);

                %Nondimensionalize from the start?
                if B1 == 0
                    VxRim = VxRim/(B2/2);
                    VyRim = VyRim/(B2/2);
                else
                    VxRim = VxRim/(B1/2);
                    VyRim = VyRim/(B1/2);
                end

                %Enclosure Blob Coordinates.
                %Enclosure is a ring centered about the lab origin.
                [x_Enc, y_Enc] = DiscretizeEnclosure(R,s); 

                %Translate beast geometric center to the chosen initial position in enclosure frame.
                xcoord = xcoord + x_o;
                ycoord = ycoord + y_o;

                %Grab the head coordinate for plotting
                x_head = xcoord(end - NRim + 1);
                y_head = ycoord(end - NRim + 1);

                [Ux_history, Uy_history, W_history, theta_history, x_cm_history, y_cm_history, separation_history, dt_history, time_history]...
                    = TimeAdvance(T, dt, xcoord, ycoord, x_Enc, y_Enc, theta_o, epsilon, VxRim, VyRim, NRim, R, a, Scale);

                speed_history = (Uy_history.^2 + Ux_history.^2).^(1/2);
                NEnc = length(y_Enc);
                
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

                %% Plot the cm trajectory
                fig = figure(1);
                ax1 = axes('Position',[0 0 1 1],'Visible','off');
                ax2 = axes('Position',[.3 .1 .6 .8]);
                axes(ax1);
                text(.025,0.6,descr);
                %back to data
                axes(ax2);

                scatter(x_cm_history(1), y_cm_history(1), '.', 'g'); %Begin at green
                hold on
                scatter(x_cm_history(end-1), y_cm_history(end-1), '.', 'r'); %End at red
                scatter(x_cm_history(2:end-2), y_cm_history(2:end-2), '.', 'k');
                daspect([1,1,1])
                scatter(x_Enc, y_Enc, '.', 'b');
                axis off
                hold off
                title1 = ['trajectory' num2str(i) num2str(j) num2str(k) num2str(m)];
                saveas(gcf,title1, 'png')
                %% Plotting velocity vs separation
                fig = figure(2);
                ax1 = axes('Position',[0 0 1 1],'Visible','off');
                ax2 = axes('Position',[.3 .1 .6 .8]);
                axes(ax1);
                text(.025,0.6,descr);
                %back to plotting data
                axes(ax2);
                hold on
                scatter(separation_history(2:end-1), speed_history(2:end-1));
                title('Swimming Speed vs. Distance from Enclosure Wall','FontSize',16,'FontWeight','bold')
                xlabel('Nondimensional Separation Distance (a/l_s)')
                ylabel('Nondimensional Swimming Speed')
                hold off
                title2 = ['speed' num2str(i) num2str(j) num2str(k) num2str(m)];
                saveas(gcf,title2, 'png')
                %% Plotting angular velocity vs separation
                fig = figure(3);
                ax1 = axes('Position',[0 0 1 1],'Visible','off');
                ax2 = axes('Position',[.3 .1 .6 .8]);
                axes(ax1);
                text(.025,0.6,descr);
                %back to plotting data
                axes(ax2);
                hold on
                scatter(separation_history(2:end-1), W_history(2:end-1));
                title('Angular velocity vs. Distance from Enclosure Wall','FontSize',16,'FontWeight','bold')
                xlabel('Nondimensional Separation Distance (a/l_s)')
                ylabel('Nondimensional Angular Velocity')
                hold off
                title3 = ['angular' num2str(i) num2str(j) num2str(k) num2str(m)];
                saveas(gcf,title3, 'png')
                %% Plotting separation vs simulation time
                fig = figure(4);
                ax1 = axes('Position',[0 0 1 1],'Visible','off');
                ax2 = axes('Position',[.3 .1 .6 .8]);
                axes(ax1);
                text(.025,0.6,descr);
                %back to plotting data
                axes(ax2);
                scatter(time_history(2:end-1), separation_history(2:end-1));
                title('Distance from Enclosure Wall vs Time','FontSize',16,'FontWeight','bold')
                xlabel('Nondimensional Time')
                ylabel('Nondimensional Separation')
                title4 = ['separation' num2str(i) num2str(j) num2str(k) num2str(m)];
                saveas(gcf,title4, 'png')
            end
        end
    end
end