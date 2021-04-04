function [] = iterPolarplot(fvalue, x, y, min, max)
    % José Reis, 2021
    %{
        fvalue  = frequency range values (array)
        x       = theta or phi values (array)
        y       = Frequency response (array)
        min     = starting position of x
        max     = ending position of x
    %}

    f = figure;
    % Retrieve figure's width
    f_w = f.Position(3);
    
    % Initialize plot's figure.
    polarplot(x, y(:,min), 'Color', 'b', 'LineWidth', 1.5)
    
    % Generate GUI Slidebar, along with its labels
    % (SB = slidebar | L_lbl = Left label | R_lbl = Right Label)
    % (M_LBL = frequency indicator label)
    sb = uicontrol( 'Parent',f,...
                    'Style','slider',...
                    'Position',[81,0,419,25],...
                    'value', min,...
                    'min', min,...
                    'max', max);

    L_lbl = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[50,0,25,20],...
                        'String',num2str(min),...
                        'BackgroundColor',f.Color);

    R_lbl = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[505,0,40,20],...
                        'String',"fs/2");
                    
    M_lbl = uicontrol(  'Parent',f,...
                        'Style','text',...
                        'Position',[-20,30,200,23],...
                        'String',join(['Frequency (Hz) = ', num2str(round(fvalue(2))), "Hz"]),...
                        'BackgroundColor',f.Color);
        
    % Nested function to update data from plot and label
    % using the SB value.
    function updateData
        polarplot(x,y(:,round(sb.Value)), 'Color', 'b', 'LineWidth', 1.5);
        set(M_lbl,'String', join(['Frequency (Hz) = ', num2str(round(fvalue(round(sb.Value)))), "Hz"]));
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

