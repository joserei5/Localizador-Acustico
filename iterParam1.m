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
    subplot(211)
    subplot(212)
    subplot(211)
    
    % _________________________.~Initial setup~.___________________________
    % sound velocity
    C = 20.05*sqrt(273.15+room.temp);
    
    % theta range
    theta = -180:180;
    
    % frequency range
    freqr = 50:2000;
    
    % distance range
    d = linspace(src.d.min,src.d.max,Nd); 
    % deltaX range
    d_x = linspace(room.mic.dX.min, room.mic.dX.max, Nx);
    % epsilon range
    e = linspace(epsilon.min, epsilon.max,Ne);

    recL = [-d_x/2; zeros(1,Nx)];
    recR = [+d_x/2; zeros(1,Nx)];

    source = [d'*cosd(theta(180+90+1)) d'*sind(theta(180+90+1))];
    d1 = sqrt((recL(1, :)-source(:, 1)).^2 + (recL(2, :)-source(:, 2)).^2);
    d2 = sqrt((recR(1, :)-source(:, 1)).^2 + (recR(2, :)-source(:, 2)).^2);
    
    % difference of distance in recL in relation to recR
    d_d = d1 - d2;
    d_t = d_d/C;
    
    % bar settings (GUI)
    minb_e = 1;
    maxb_e = Ne;
    minb_th = 1;
    maxb_th = 360+1;
    minb_f = 1;
    maxb_f = length(freqr);
    
    
    % ____________________________.~Algorithm~.____________________________
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

    
    % ______________________________.~Curve~.______________________________
    axf = findall(f,'type','axes');
    axf(2).PositionConstraint = 'OuterPosition';
    axf(2).Position(4) = axf(2).Position(4)+axf(2).Position(2)*0.3;
    axf(2).Position(2) = axf(2).Position(2)-axf(2).Position(2)*0.3;
    axf(1).PositionConstraint = 'OuterPosition';
    axf(1).Position(4) = axf(1).Position(4)-axf(1).Position(2)*0.8;
    axf(1).Position(2) = axf(1).Position(2)-axf(1).Position(2)*0.8;
    axf(1).Position(1) = axf(1).Position(1)+axf(1).Position(1)*5;
    axf(1).PositionConstraint = 'InnerPosition';
    axf(1).Position(3) = axf(1).Position(3)-axf(1).Position(3)*0.75;
    
    % visualize curve f(delta_x, distance, error)
    pcolor(d_x, d, abs(azimuth_e));
    tstr =  "Frequency = 500Hz | Theta = 90º | " +...
            "Epsilon = " + num2str(e(1)*1e6) + "uS";
    title(tstr)
    shading interp
    xlabel('microphone spacing (m)')
    ylabel('source distance (m)')
    zlabel('azimuth error (º)')
    zlim([0 45])
    caxis([0 45])
    colorbar
    
    % near field curve
    hold on;
    mfreq = 700;
    LAMBD = C/mfreq;
    R = (2*d.^2)./LAMBD;
    
    plot(d, R, 'r')
    hold off;
    
    % _________________________.~Head Reference~.__________________________
    axes(axf(1))
    plot(cos(0:0.05:2*pi),sin(0:0.05:2*pi),'k');
    hold on
    plot([0 cosd(90)], [0 sind(90)], 'r');
    hold off
    set(gca,'XTick',[], 'YTick', []);
    set(gca,'color','none')
    axis off
    xlim([-1 1]);ylim([-1 1]);
    disableDefaultInteractivity(axf(1))
    
    % Generate GUI Slidebar, along with its labels
    % sb_e = epsilon slidebar       | sb_th = theta slidebar
    % sb_f = frequency slidebar     
    % L_lbl_e = epsilon Left label  | R_lbl_e = epsilon Right Label
    % L_lbl_th = theta Left label   | R_lbl_th = theta Right Label
    % L_lbl_f = freq Left label     | R_lbl_f = freq Right Label
    
    sb_e = uicontrol( 'Parent',f,...
                    'Style','slider',...
                    'Position',[61,25,319,25],...
                    'value', minb_e,...
                    'min', minb_e,...
                    'max', maxb_e);
    
    sb_th = uicontrol( 'Parent',f,...
                    'Style','slider',...
                    'Position',[61,50,319,25],...
                    'value', 180+90+1,...
                    'min', minb_th,...
                    'max', maxb_th);
    
    sb_f = uicontrol( 'Parent',f,...
                    'Style','slider',...
                    'Position',[61,75,319,25],...
                    'value', 500,...
                    'min', minb_f,...
                    'max', maxb_f);

    L_lbl_e = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[5,25,45,20],...
                        'String',join([num2str(e(1)*1e6), " uS"]),...
                        'BackgroundColor',f.Color);

    L_lbl_th = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[10,50,45,20],...
                        'String',join([num2str(theta(1)), "º"]),...
                        'BackgroundColor',f.Color);

    L_lbl_f = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[10,75,45,20],...
                        'String',join([num2str(freqr(1)), "Hz"]),...
                        'BackgroundColor',f.Color);

    R_lbl_e = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[385,25,45,20],...
                        'String',join([num2str(e(Ne)*1e6), " uS"]));

    R_lbl_th = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[380,50,45,20],...
                        'String',join([num2str(theta(end)), "º"]));

    R_lbl_f = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[380,75,45,20],...
                        'String',join([num2str(freqr(end)*1e-3), "kHz"]));


    % Nested function to update data from plot and label
    % using the sb_e value.
    function updateData
        bIndex_e = round(sb_e.Value);
        actual_e = e(bIndex_e)*1e6;
        bIndex_th = round(sb_th.Value);
        actual_th = theta(bIndex_th);
        bIndex_f = round(sb_f.Value);
        actual_f = freqr(bIndex_f);
        
        source = [d'*cosd(actual_th) d'*sind(actual_th)];
        d1 = sqrt((recL(1, :)-source(:, 1)).^2 + (recL(2, :)-source(:, 2)).^2);
        d2 = sqrt((recR(1, :)-source(:, 1)).^2 + (recR(2, :)-source(:, 2)).^2);

        % difference of distance in recL in relation to recR
        d_d = d1 - d2;
        d_t = d_d/C;
        
        arg = C*(d_t+0)./d_x;
        arg(arg<-1) = -1;
        arg(arg>1) = 1;
        realAz = acosd(arg);

        arg = C*(d_t+e(bIndex_e))./d_x;
        arg(arg<-1) = -1;
        arg(arg>1) = 1;
        azimuth = acosd(arg);
        azimuth_e = realAz - azimuth;

        % ____________________________.~Curve~.____________________________
        % visualize curve f(delta_x, distance, error)
        axes(axf(2))
        pcolor(d_x, d, abs(azimuth_e));
        shading interp
        tstr =  "Frequency = " + num2str(actual_f) + "Hz | " +...
                "Theta = " + num2str(actual_th) + " º | " +...
                "Epsilon = " + num2str(actual_e) + "uS";
        title(tstr)
        xlabel('microphone spacing (m)')
        ylabel('source distance (m)')
        zlabel('azimuth error (º)')
        zlim([0 45])
        caxis([0 45])
        colorbar
        
        % near field curve
        hold on;
        LAMBD = C/actual_f;
        R = (2*d.^2)./LAMBD;
        
        plot(d, R, 'r')
        hold off;
        
        % _________________________.~Head Reference~.__________________________
        axes(axf(1))
        plot(cos(0:0.05:2*pi),sin(0:0.05:2*pi),'k');
        hold on
        plot([0 0+cosd(actual_th)], [0 sind(actual_th)], 'r');
        hold off
        set(gca,'XTick',[], 'YTick', []);
        set(gca,'color','none')
        axis off
        
    end

    % Nested function to update the GUI resolution using
    % the figure's position values.
    function updateGUI
        % Calculate difference of width
        diff_w = f.Position(3) - f_w;
        
        % Estimate SB and Label's new width
        sb_w = 319 + diff_w;
        rl_w = 385 + diff_w;
        
        % Truncate min width to 0
        if(sb_w < 0) sb_w = 0; end
        if(rl_w < 0) rl_w = 0; end
        
        % Update GUI resolution
        set(sb_e, 'Position', [61, 25, sb_w, 25]);
        set(sb_th, 'Position', [61, 50, sb_w, 25]);
        set(sb_f, 'Position', [61, 75, sb_w, 25]);
        set(R_lbl_e, 'Position', [rl_w, 25, 45, 20]);
        set(R_lbl_th, 'Position', [rl_w-5, 50, 45, 20]);
        set(R_lbl_f, 'Position', [rl_w-5, 75, 45, 20]);
    end
    
    % React to changes in the SB or figure
    addlistener(sb_e, 'ContinuousValueChange', @(es,ed) updateData);
    addlistener(sb_th, 'ContinuousValueChange', @(es,ed) updateData);
    addlistener(sb_f, 'ContinuousValueChange', @(es,ed) updateData);
    f.SizeChangedFcn = @(es,ed)updateGUI;
end

