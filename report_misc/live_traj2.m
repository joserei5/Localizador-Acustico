clear all;clc;close all;
addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generic
addpath ../structures

%%
global trjF;
trjF = figure();
trjF.ToolBar = 'none';
trjF.MenuBar = 'none';
trjF.Resize = 'off';
set(trjF, 'WindowButtonDownFcn', @PrintValues);
% set(trjF, 'WindowButtonMotionFcn', @Crosshair);

global px;
global py;
global trjP;
trjP = plot(1,1);
trjP.XDataSource = 'px';
trjP.YDataSource = 'py';
xlim([0 10]);
ylim([0 10]);

global chp;
chp = plot(1,1);

function Crosshair(~,~)
    global trjP;
    hold on;
    P = get (gca, 'CurrentPoint');
    x = P(1,1);
    y = P(1,2);
    global chp;
    chp = plot(x,y);
    xline(x);
    yline(y);
    hold off;
end

function PrintValues(object,event)
    global trjP;
    global px;
    global py;
    P = get (gca, 'CurrentPoint');
    x=P(1,1);
    y=P(1,2);
%     if (x>1 && x<2)
%         if (y>1 && y<2)
            px=[px x];
            py=[py y];
            refreshdata;
            drawnow;
            fprintf('(%.2f , %.2f)\n',x,y)
%         end
%     end

    disp(get(gcf,'KeyPressFcn'));
end