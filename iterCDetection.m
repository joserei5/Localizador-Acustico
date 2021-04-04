function [] = iterCDetection(room, rec, distance, points)
    % Iterative Circle Detection
    % José Reis, 2021
    %{
        room        = room structure
        rec         = receiver structure
        distance    = trajectory distance to receiver
        points      = number of points to detect
    %}

    f = figure;
    % Retrieve figure's width
    f_w = f.Position(3);
    
    offset = 90;
    [rp, th_azimuth] = speaker_circle2D(room.size, ...
                                        room.rec.loc,...
                                        distance,...
                                        points, offset);
    th_azimuth = -th_azimuth + offset + 360;

    minb = 1;
    maxb = length(rp);
    
    [sndfile, fs] = audioread('cportugal.wav');
    rangemax = round(rec.mic.dist + 1e-3*fs);
    
    C = 20.05*sqrt(273.15+room.temp);
    
    % Initialize plot's figure.
    px = room.size(minb,1);
    py = room.size(minb,2);
    pz = room.size(minb,3);
    X = [0 px px 0; 0 px px 0; ...
         0 0 px px; 0 0 px px;
         0 0 0 0; px px px px];
    Y = [0 0 py py; 0 0 py py; ...
         0 0 0 0; py py py py;
         0 0 py py; 0 0 py py];
    Z = [0 0 0 0; pz pz pz pz;
         0 pz pz 0; 0 pz pz 0; 
         0 pz pz 0; 0 pz pz 0];
    pcol = [0.5 0.5 0.7];
    
    for i=1:size(X,1)
        fill3(X(i,:), Y(i,:), Z(i,:), pcol+[0.05*i 0.05*i 0], 'facealpha', 0.2)
        hold on
    end
    axis([0 px 0 py 0 pz]);
    daspect([1 1 1])
    plot3(rp(minb,1), rp(minb,2), rp(minb,3), 'b*', 'LineWidth', 2)
    
    ox = room.rec.loc(1);
    oy = room.rec.loc(2);
    oz = room.rec.loc(3);
    plot3(ox, oy, oz, 'k+')
    
    x = cos(linspace(0,2*pi,100))*distance + ox;
    y = sin(linspace(0,2*pi,100))*distance + oy;
    z = ones(1,100) * oz;
    plot3(x,y,z,'k-')
    
    % draw predicted and theoretical line
    point1 = [0 rec.mic.dist 0];
    point2 = [0 -rec.mic.dist 0];
    point1 = rotate2D(point1, rec.yaw);
    point2 = rotate2D(point2, rec.yaw);
    lx = [point1(1) point2(1)]+ox;
    ly = [point1(2) point2(2)]+oy;
    lz = [point1(3) point2(3)]+oz;
    plot3(lx, ly, lz, 'b--')

    [yL, yR] = sim_stereo(  rec.struct,... 
                            room.size,...
                            room.rec.loc,...
                            [rec.yaw, rec.pitch],...
                            rp(minb,:),...
                            sndfile, fs,...
                            1,...
                            rec.mic.dmf);

    lr_corr = xcorr(yL,yR, rangemax);
    [value, index] = max(abs(lr_corr));
    delay_index = (index-1) - rangemax;
    delay_t = delay_index*1/fs;
    arg = C*delay_t/rec.mic.dist;
    if arg > 1
        arg = 1;
    elseif arg < -1
        arg = -1;
    end
    azimuth = acosd(arg);
    e_azimuth = th_azimuth(1) - azimuth;

    hold on;
    pL = sqrt((ox-rp(minb,1))^2 + (oy-rp(minb,2))^2);
    p0 = room.rec.loc;
    DAz = 90 - azimuth;
    aX = p0(1)+pL*(cosd(DAz));
    aY = p0(2)-pL*(sind(DAz));
    p1 = rp(minb,:);
    p2 = [aX aY oz];

    LX1 = [p0(1) p1(1)];
    LY1 = [p0(2) p1(2)];
    LZ1 = [p0(3) p1(3)];

    LX2 = [p0(1) p2(1)];
    LY2 = [p0(2) p2(2)];
    LZ2 = [p0(3) p2(3)];

    plot3(LX1, LY1, LZ1, 'k');
    plot3(LX2, LY2, LZ2, 'r');
    
    view(2)
    
    hold off
    
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
                        'String',num2str(minb),...
                        'BackgroundColor',f.Color);

    R_lbl = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[505,0,40,20],...
                        'String',"180º");
    
    Title_lbl = uicontrol(  'Parent',f,...
                            'Style','text',...
                            'Position',[20,90,100,23],...
                            'String',"YAW:");
                    
    M_lbl = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[0,70,140,23],...
                        'String',join(['Predicted = ', num2str(azimuth), "º"]),...
                        'BackgroundColor',f.Color);

    M_lbl2 = uicontrol( 'Parent',f,...
                        'Style','text',...
                        'Position',[0,50,140,23],...
                        'String',join(['Theoretical = ', num2str(th_azimuth(1)), "º"]),...
                        'BackgroundColor',f.Color);

    M_lbl3 = uicontrol( 'Parent',f,...
                        'Style','text',...
                        'Position',[0,30,140,23],...
                        'String',join(['Error = ', num2str(e_azimuth), "º"]),...
                        'BackgroundColor',f.Color);
        
    % Nested function to update data from plot and label
    % using the SB value.
    function updateData
        bIndex = round(sb.Value);
        
        px = room.size(1);
        py = room.size(2);
        pz = room.size(3);
        X = [0 px px 0; 0 px px 0; ...
             0 0 px px; 0 0 px px;
             0 0 0 0; px px px px];
        Y = [0 0 py py; 0 0 py py; ...
             0 0 0 0; py py py py;
             0 0 py py; 0 0 py py];
        Z = [0 0 0 0; pz pz pz pz;
             0 pz pz 0; 0 pz pz 0; 
             0 pz pz 0; 0 pz pz 0];
        pcol = [0.5 0.5 0.7];

        for i=1:size(X,1)
            fill3(X(i,:), Y(i,:), Z(i,:), pcol+[0.05*i 0.05*i 0], 'facealpha', 0.2)
            hold on
        end
        axis([0 px 0 py 0 pz]);
        daspect([1 1 1])

        plot3(rp(bIndex,1), rp(bIndex,2), rp(bIndex,3), 'b*', 'LineWidth', 2)
        plot3(ox, oy, oz, 'k+')
        
        x = cos(linspace(0,2*pi,100))*distance + ox;
        y = sin(linspace(0,2*pi,100))*distance + oy;
        z = ones(1,100) * oz;
        plot3(x,y,z,'k-')
        
        % draw predicted and theoretical line
        point1 = [0 rec.mic.dist 0];
        point2 = [0 -rec.mic.dist 0];
        point1 = rotate2D(point1, rec.yaw);
        point2 = rotate2D(point2, rec.yaw);
        lx = [point1(1) point2(1)]+ox;
        ly = [point1(2) point2(2)]+oy;
        lz = [point1(3) point2(3)]+oz;
        plot3(lx, ly, lz, 'b--')

        [yL, yR] = sim_stereo(  rec.struct,... 
                                room.size,...
                                room.rec.loc,...
                                [rec.yaw, rec.pitch],...
                                rp(bIndex,:),...
                                sndfile, fs,...
                                1,...
                                rec.mic.dmf);

        lr_corr = xcorr(yL,yR, rangemax);
        [value, index] = max(abs(lr_corr));
        delay_index = (index-1) - rangemax;
        delay_t = delay_index*1/fs;
        arg = C*delay_t/rec.mic.dist;
        if arg > 1
            arg = 1;
        elseif arg < -1
            arg = -1;
        end
        azimuth = acosd(arg);
        e_azimuth = th_azimuth(bIndex) - azimuth;

        hold on;
        pL = sqrt((ox-rp(bIndex,1))^2 + (oy-rp(bIndex,2))^2);
        p0 = room.rec.loc;
        DAz = 90 - azimuth;
        aX = p0(1)+pL*(cosd(DAz));
        aY = p0(2)-pL*(sind(DAz));
        p1 = rp(bIndex,:);
        p2 = [aX aY oz];

        LX1 = [p0(1) p1(1)];
        LY1 = [p0(2) p1(2)];
        LZ1 = [p0(3) p1(3)];

        LX2 = [p0(1) p2(1)];
        LY2 = [p0(2) p2(2)];
        LZ2 = [p0(3) p2(3)];

        plot3(LX1, LY1, LZ1, 'k');
        plot3(LX2, LY2, LZ2, 'r');
        
%         angle = 90 - acosd((rp(bIndex,1) - ox)/distance);
        set(M_lbl,'String', join(['Predicted = ', num2str(azimuth), "º"]));
        set(M_lbl2,'String', join(['Theoretical = ', num2str(th_azimuth(bIndex)), "º"]));
        set(M_lbl3,'String', join(['Error = ', num2str(e_azimuth), "º"]));
        hold off
        view(2)
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

