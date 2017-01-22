% Plots the trajectory for a boundary element squirmer in an enclosed domain 

tic


%These are scaled to become the default simulation time T and timestep
%dt_o, although both will change depending on how dt is rescaled during the
%simulation.
Time = 35;
Increment = 0.5;

%When h ~=0 an if statement replaces the circular enclosure with a vertical
%wall

theta_o_array = [pi/10 -pi/10 pi/7 -pi/7 pi/4 -pi/4 pi/3 -pi/3 0 pi];
r_o_array = [0.5 0.7 0.85]; %Ratio r_o/R
a_array = [1e-3 1e-1 1e1]; % radius of the disk nondimensionalized by the Saffman length
B2_array = [-1 1];

B1 = 0;

%that is placed symmetrically at a distance R to the
%right of the origin

FinishedLoopCount = 0;

Steps = floor(Time/Increment);
x_cm_trajectories = zeros([Steps+1,length(theta_o_array)]);
y_cm_trajectories = zeros([Steps+1,length(theta_o_array)]);



for i = 1:length(a_array) %Cycle through initial orientation
            R = 10; %Place enclosure R*a distance to the right of the origin
            h = 15; % Ratio h/a
            
            a = a_array(i);  %%% radius of the disk nondimensionalized by the Saffman length
            R = R*a;  %Rescale in terms of a
            h = h*a;
            
            s = 0.12 * a;          %%% radial spacing between neighboring blobs
            epsilon = s/8;        %%% radius of blobs

            %%%Needs to be implemented into solve_U_enclosure
            Scale = 10/a;         %%% This scales some terms in the matrix equation. Possibly helps with invertibility.

            %This rescaling makes the trajectories similar in distance
            %when compared to the enclosure size. (This way I can cycle
            %through many beast radii and leave [Time] fixed.)
            T = Time*a;         
            dt = Increment*a;  
            
    for j = 1:length(B2_array) %Cycle through initial separation
            B2 = B2_array(j);
            
        for k = 1:length(r_o_array) %Cycle reduced radius
            r_o = r_o_array(k)*R;   %Place beast radial distance r_o*R from the origin
                            
            for m = 1:length(theta_o_array) %Cycle initial orientation
                close all
                
                theta_o = theta_o_array(m);   %%% Beast intial orientation (head direction in lab frame)
                                

                

                %Initial Conditions

                phi_o = 0;       %%% Angle coordinate of beast cm from origin
                
                x_o = r_o*cos(phi_o); %%% Beast CM initial x position as seen in lab frame.
                y_o = r_o*sin(phi_o); %%% Beast CM Initial y position in lab frame.




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
                LoopPercentage = FinishedLoopCount/(length(a_array)*length(B2_array)*length(r_o_array)*length(theta_o_array))*100
                
                x_cm_trajectories(:,m) = x_cm_history;
                y_cm_trajectories(:,m) = y_cm_history;

                toc
            end
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
                    str_B2;
                    str_a;
                    str_r_o;
                    str_R};

                %% Plot the cm trajectory
                fig = figure(1);
                ax1 = axes('Position',[0 0 1 1],'Visible','off');
                ax2 = axes('Position',[.3 .1 .6 .8]);
                axes(ax1);
                text(.1,0.6,descr);
                %back to data
                axes(ax2);

                
                hold on
                
                for n = 1:m
                    %Plot CM
                    scatter(x_cm_trajectories(:,n), y_cm_trajectories(:,n), '.', 'k');
                    scatter(x_cm_trajectories(1,n), y_cm_trajectories(1,n),'o','g'); %Begin at green
                    scatter(x_cm_trajectories(end,n), y_cm_trajectories(end,n),'o','r'); %End at red
                    
                    %Body. Inner then rim. At two points in time (1 and end)
                    %scatter(x_history(end,1:Nbeast-NRim),y_history(end,1:Nbeast-NRim),'.','b')
                    %scatter(x_history(end,Nbeast-NRim+1:end),y_history(end,Nbeast-NRim+1:end),'.','r')
                    %scatter(x_history(1,1:Nbeast-NRim),y_history(1,1:Nbeast-NRim),'.','b')
                    %scatter(x_history(1,Nbeast-NRim+1:end),y_history(1,Nbeast-NRim+1:end),'.','r')

                    %Head
                    %scatter(x_history(end,Nbeast - NRim + 1),y_history(end,Nbeast - NRim + 1),400,'.','m')
                    %scatter(x_history(1,Nbeast - NRim + 1),y_history(1,Nbeast - NRim + 1),400,'.','m')

                    
                end
                

                
                
                daspect([1,1,1])
                scatter(x_Enc, y_Enc, '.','b');
                axis off
                hold off
                title1 = ['trajectory' num2str(i) num2str(j) num2str(k) num2str(m)];
                saveas(gcf,title1, 'png')
                saveas(gcf,title1, 'eps')
                saveas(gcf,title1, 'fig')
        end
    end
end