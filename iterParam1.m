function [] = iterParam1(room, src, epsilon, Nx, Nd, Ne)
    % Iterative Parameter simulator
    % José Reis, 2021
    %{
        room        = room structure
        src         = source structure
        epsilon     = epsilon range
        points      = number of points to simulate
    %}

    f = figure;
    % Retrieve figure's width
    f_w = f.Position(3);
    
    % _________________________.~Initial setup~.___________________________
    % sound velocity
    C = 20.05*sqrt(273.15+room.temp);
    
    % update reference (90º = 0º)
    src.theta = src.theta+90;
    
    % distance range
    d = linspace(src.d.min,src.d.max,Nd); 
    % deltaX range
    d_x = linspace(room.mic.dX.min, room.mic.dX.max, Nx);
    % epsilon range
    e = linspace(epsilon.min, epsilon.max,Ne);
    
%     % source point
%     source = [d*cosd(src.theta) d*sind(src.theta)];
%     % recL point
%     recL = [-d_x/2 zeros(Nx, 1)];
%     % recR point
%     recR = [+d_x/2 zeros(Nx, 1)];
%     
%     % calculate distance to recL and recR
%     d1 = sqrt((recL(:,1)-source(1)).^2 + (0-source(2)).^2);
%     d2 = sqrt((recR(:,1)-source(1)).^2 + (0-source(2)).^2);

    d_d = zeros(Nd,Nx);
    
%     for nd=1:Nd
%         source = [d(nd)*cosd(src.theta) d(nd)*sind(src.theta)];
%         for nx=1:Nx
%             recL = [-d_x(nx)/2 0];
%             recR = [+d_x(nx)/2 0];
%             d1 = sqrt((recL(1)-source(1)).^2 + (recL(2)-source(2)).^2);
%             d2 = sqrt((recR(1)-source(1)).^2 + (recR(2)-source(2)).^2);
%             d_d(nd, nx) = d1 - d2;
%         end
%     end

    recL = [-d_x/2; zeros(1,Nx)];
    recR = [+d_x/2; zeros(1,Nx)];
%     for nd=1:Nd
%         source = [d(nd)*cosd(src.theta) d(nd)*sind(src.theta)];
%         d1 = sqrt((recL(1, :)-source(1)).^2 + (recL(2,:)-source(2)).^2);
%         d2 = sqrt((recR(1, :)-source(1)).^2 + (recR(2,:)-source(2)).^2);
%         d_d(nd,:) = d1-d2;
%     end

    source = [d'*cosd(src.theta) d'*sind(src.theta)];
    d1 = sqrt((recL(1, :)-source(:, 1)).^2 + (recL(2, :)-source(:, 2)).^2);
    d2 = sqrt((recR(1, :)-source(:, 1)).^2 + (recR(2, :)-source(:, 2)).^2);
    
    % difference of distance in recL in relation to recR
    d_d = d1 - d2;
    d_t = d_d/C;
    
    % prepare variables for 3D curve
%     [d_x, d] = meshgrid(d_x, d);
    
    % bar settings (GUI)
    minb = 1;
    maxb = Ne;
    
    % original azimuth
%     SOL = ones(Nx,1)*src.theta;
    
    
    % ____________________________.~Algorithm~.____________________________
%     arg = C*(d_t+0)./d_x;
%     arg(arg<-1) = -1;
%     arg(arg>1) = 1;
%     realAz = acosd(arg);
%     
%     arg = C*(d_t+e(1))./d_x;
%     arg(arg<-1) = -1;
%     arg(arg>1) = 1;
%     azimuth = acosd(arg);
% %     azimuth_e = SOL - azimuth;
%     azimuth_e = realAz - azimuth;

    realAz = zeros(Nd,Nx);      % ref for \epsilon=0
    azimuth = zeros(Nd,Nx);     % actual values
    azimuth_e = zeros(Nd,Nx);   % azimuth error
    
    arg = C*(d_t+0)./d_x;
    arg(arg<-1) = -1;
    arg(arg>1) = 1;
    realAz = acosd(arg);

    arg = C*(d_t+e(1))./d_x;
    arg(arg<-1) = -1;
    arg(arg>1) = 1;
    azimuth = acosd(arg);
    azimuth_e = realAz - azimuth;
    
%     for nd=1:Nd
%         for nx=1:Nx
%             arg = C*(d_t(nd,nx)+0)./d_x(nx);
%             arg(arg<-1) = -1;
%             arg(arg>1) = 1;
%             realAz(nd,nx) = acosd(arg);
%             
%             arg = C*(d_t(nd,nx)+e(1))./d_x(nx);
%             arg(arg<-1) = -1;
%             arg(arg>1) = 1;
%             azimuth(nd,nx) = acosd(arg);
%             azimuth_e(nd,nx) = realAz(nd,nx) - azimuth(nd,nx);
%         end
%     end
    
    % ______________________________.~Curve~.______________________________
    % visualize curve f(delta_x, distance, error)
    pcolor(d_x, d, abs(azimuth_e));
    tstr = "Epsilon = " + num2str(e(1)*1e6) + "uS";
    title(tstr)
    shading interp
    xlabel('microphone spacing (m)')
    ylabel('source distance (m)')
    zlabel('azimuth error (º)')
    zlim([0 45])
    caxis([0 45])
    colorbar
    
    % Generate GUI Slidebar, along with its labels
    % (SB = slidebar | L_lbl = Left label | R_lbl = Right Label)
    % (M_LBL = frequency indicator label)
    sb = uicontrol( 'Parent',f,...
                    'Style','slider',...
                    'Position',[81,0,419,25],...
                    'value', minb,...
                    'min', minb,...
                    'max', maxb);

    L_lbl = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[50,0,25,20],...
                        'String',join([num2str(e(1)*1e6), " uS"]),...
                        'BackgroundColor',f.Color);

    R_lbl = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[505,0,40,20],...
                        'String',join([num2str(e(Nx)*1e6), " uS"]));
                    
        
    % Nested function to update data from plot and label
    % using the SB value.
    function updateData
        bIndex = round(sb.Value);
        actual_e = e(bIndex)*1e6;
        
        % __________________________.~Algorithm~.__________________________
%         arg = C*(d_t+e(bIndex))./d_x;
%         arg(arg<-1) = -1;
%         arg(arg>1) = 1;
%         azimuth = acosd(arg);
% %         azimuth_e = SOL - azimuth;
%         azimuth_e = realAz - azimuth;

%         for nd=1:Nd
%             for nx=1:Nx
%                 arg = C*(d_t(nd,nx)+e(bIndex))./d_x(nx);
%                 arg(arg<-1) = -1;
%                 arg(arg>1) = 1;
%                 azimuth = acosd(arg);
%                 azimuth_e(nd,nx) = realAz(nd,nx) - azimuth;
%             end
%         end

        arg = C*(d_t+0)./d_x;
        arg(arg<-1) = -1;
        arg(arg>1) = 1;
        realAz = acosd(arg);

        arg = C*(d_t+e(bIndex))./d_x;
        arg(arg<-1) = -1;
        arg(arg>1) = 1;
        azimuth = acosd(arg);
        azimuth_e = realAz - azimuth;

        % ____________________________.~Curve~.____________________________
        % visualize curve f(delta_x, distance, error)
        pcolor(d_x, d, abs(azimuth_e));
        shading interp
        tstr = "Epsilon = " + num2str(actual_e) + "uS";
        title(tstr)
        xlabel('microphone spacing (m)')
        ylabel('source distance (m)')
        zlabel('azimuth error (º)')
        zlim([0 45])
        caxis([0 45])
        colorbar
    end

    % Nested function to update the GUI resolution using
    % the figure's position values.
    function updateGUI
        % Calculate difference of width
        diff_w = f.Position(3) - f_w;
        
        % Estimate SB and Label's new width
        sb_w = 419 + diff_w;
        rl_w = 505 + diff_w;
        
        % Truncate min width to 0
        if(sb_w < 0) sb_w = 0; end
        if(rl_w < 0) rl_w = 0; end
        
        % Update GUI resolution
        set(sb, 'Position', [81, 0, sb_w, 25]);
        set(R_lbl, 'Position', [rl_w, 0, 40, 20]);
    end
    
    % React to changes in the SB or figure
    addlistener(sb, 'ContinuousValueChange', @(es,ed) updateData);
    f.SizeChangedFcn = @(es,ed)updateGUI;
end

