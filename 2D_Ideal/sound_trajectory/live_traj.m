clear all;clc;close all;
addpath ../../simLocUS
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../structures

%% Variables
global room
room.xy = [5 7];
room.rc = [0 0 0 0];

global rec;
load('recstruct.mat');
rec.type = 0;
rec.loc = [1 3.5 1.7];
rec.mic.dmf = 1;
rec.mic.dist = 0.1*rec.mic.dmf*2;
rec.th = 0;
rec.phi = 0;

spk.loc=rec.loc;

global sDrawing;
global dDrawing; dDrawing=0;
sDrawing=0;

global tX;
global tY;

global WAV;
WAV = 'mosquito.wav';
global v;
v = 0.5;


%% Initialize figure
global trjF;
trjF = figure();
trjF.ToolBar = 'none';
trjF.MenuBar = 'none';
trjF.Resize = 'off';
set(gcf, 'KeyPressFcn', @gui);
% set(gcf, 'WindowButtonMotionFcn', @pencil);
% set(gcf, 'WindowButtonDownFcn', @actuator);

% hold on;
% global tPlot;
% tPlot = plot(tX, tY, 'k');
% hold off;

%% Receiver
% hold on;
% rec_offset=rec.mic.pos(1,2);
% plot(rec.loc(1), rec.loc(2), 'r+')
% plot(rec.loc(1), rec.loc(2)+rec_offset, 'r>')
% plot(rec.loc(1), rec.loc(2)-rec_offset, 'r>')
% xlim([0 room.xy(1)])
% ylim([0 room.xy(2)])
% hold off;

%% Detection Lines
% hold on;
% % detection line
% det_src_x = spk.loc(1,1);
% det_src_y = spk.loc(1,2);
% det_line_x = [rec.loc(1) det_src_x];
% det_line_y = [rec.loc(2) det_src_y];
% det_plt= plot(det_line_x, det_line_y, 'r');
% det_plt.XDataSource = 'det_line_x';
% det_plt.YDataSource = 'det_line_y';
% 
% % theoretical azimuth line
% taz_src_x = spk.loc(1,1);
% taz_src_y = spk.loc(1,2);
% taz_line_x = [rec.loc(1) taz_src_x];
% taz_line_y = [rec.loc(2) taz_src_y];
% taz_plt= plot(taz_line_x, taz_line_y, 'k');
% taz_plt.XDataSource = 'taz_line_x';
% taz_plt.YDataSource = 'taz_line_y';
% hold off;

%% Functions
function pencil (object, eventdata)
    P = get (gca, 'CurrentPoint');
    Px=P(1,1);
    Py=P(1,2);
    global room;
    global sDrawing;
    
    global trjF;
    global tX;
    global tY;
    global ttText;
    
    if sDrawing == 1
        hold on;
        if (Px>=0) && (Px<=room.xy(1))
            if (Py>=0) && (Py<=room.xy(2))
                tX = [tX Px]; tY = [tY Py];
                tPlot = plot(tX, tY, 'k');
            end
        end
    end
    
end

function actuator (object, eventdata)
    global sDrawing;
    global dDrawing;
    global trjF;
    global tX;
    global tY;
    
    set(gcf, 'KeyPressFcn', []);
    
    if sDrawing == 1
        sDrawing = 0;
        plot(tX(1), tY(1), 'k+')
        plot(tX(end), tY(end), 'k|')
        hold off;
        set(gcf, 'WindowButtonMotionFcn', []);
        set(gcf, 'WindowButtonDownFcn', []);
        dDrawing = 1;
        savetraj(tX,tY);
    elseif sDrawing == 0
        sDrawing = 1;
    end
end

function savetraj (X,Y)
    TRAJ.X = X;
    TRAJ.Y = Y;
    path = '/home/jreis/Localizador-Acustico/structures/trajectories/';
    fname = join(['custom',datestr(now,'ddmmyy_HHMM')]);
    fpath = join([path,fname]);
    save(fpath,'TRAJ');
end

function gui(~,event)
    if(event.Key == 's')
        drawContainer();
        set(gcf, 'WindowButtonMotionFcn', @pencil);
        set(gcf, 'WindowButtonDownFcn', @actuator);
    end
end

function drawContainer ()
    global trjF;
    figure(trjF);
    
    global room;
    rectangle('Position', [0 0 room.xy(1) room.xy(2)])
    
    global rec;
    hold on;
    rec_offset=rec.mic.pos(1,2);
    plot(rec.loc(1), rec.loc(2), 'r+')
    plot(rec.loc(1), rec.loc(2)+rec_offset, 'r>')
    plot(rec.loc(1), rec.loc(2)-rec_offset, 'r>')
    hold off;

    axis equal
    xlim([0 room.xy(1)])
    ylim([0 room.xy(2)])
    
    hold on;
    global tPlot;
    global tX;
    global tY;
    tPlot = plot(tX, tY, 'k');
    hold off;
end

function buildSound ()
    % read WAV audio
    [s, fs]=    audioread(WAV);
    Ns=         length(s);
    t=          (1:Ns)*1/fs;
    
    % frequency limits
    fl= 100;   % lower frequency bound
    fh= 20e3;  % upper frequency bound (<=80 kHz due to @KemoL10_TF)
    
    % velocity
    global v;

    % first and final point position (y-axis)
    p_x=    3;
    p_yi=   0.1;
    p_yf=   dSize(2)-0.1;

    % trajectory jump
    df=     v/fs;

    % rebuild trajectory
    p_y=    p_yi:df:p_yf;
    N=      length(p_y);

    % initialize stereo channels
    y=  zeros(2,N);

    
end