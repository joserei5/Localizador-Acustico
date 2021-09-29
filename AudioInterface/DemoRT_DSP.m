fclose('all');
close all;
clear force all;
clc;
addpath ../functions

%% Get drivers
info = audiodevinfo;
driversIN = cell(length(info.input),1);
for i=1:length(info.input)
    driversIN{i,1} = join([string(info.input(i).ID), ' - ', info.input(i).Name]);
end
driversIN

driversOUT = cell(length(info.output),1);
for i=1:length(info.output)
    driversOUT{i,1} = join([string(info.output(i).ID), ' - ', info.output(i).Name]);
end
driversOUT

%% Prepare audio devices
qc.fs = 48e3;   % sampling frequency
qc.b = 24;      % audio bits
qc.dev.r = 1;   % 1-2 (QUAD-CAPTURE) (Windows DirectSound) -- input
qc.dev.p = 9;   % 1-2 (QUAD-CAPTURE) (Windows DirectSound) -- output

BLK_t = 100*1e-3;   % Block size = 100 ms
BLK_N = BLK_t*qc.fs;% Block size = samples

deviceReader = audioDeviceReader;
deviceReader.Device = '1-2 (QUAD-CAPTURE)';
deviceReader.NumChannels = 2;
deviceReader.SampleRate = qc.fs;
deviceReader.BitDepth = join([num2str(qc.b),'-bit integer']);
deviceReader.SamplesPerFrame = BLK_N;

setup(deviceReader)

% player -- drone noise
[drone, dfs] = audioread('../soundfiles/generic/drone1.mp3');
droneObj = audioplayer(drone, dfs, qc.b, qc.dev.p);

% player -- frequency beep
t = 0:1/qc.fs:BLK_t;
beepf = 1200;
fbeep = 0.25*sin(2*pi*beepf*t);
beepObj = audioplayer(fbeep, qc.fs, qc.b, qc.dev.p);

%% Algorithm Variables
MED_t = 1*1e-3;                 % maximum estimated delay = 1 ms
MED_N = MED_t * qc.fs;          % maximum estimated delay (in samples)

ROOM.T = 20;                    % room temperature (ºC)
C = 20.05*sqrt(273.15 + ROOM.T);% sound velocity (m/s)
REC.d = 29.5e-2;                % receiver distance = 29.5 cm

CH.fs = qc.fs;
CR = round(MED_N); % correlation range=receiver distance (1 ms estimation)

AOA = [];

%% Figures
%"shed" + little space
figure;
shed = rectangle('Position',[0 0 1 1]);
shed.FaceColor = [.777 .777 .777];
hold on;

%receiver
r.x=0.7;
r.y=0.5;
plot(r.x, r.y, 'r+');

%detection line
d = sqrt(1.5^2+1.5^2);
d = 2*d;
xx=r.x;
yy=r.y;
SP1=plot(xx, yy, 'r');
SP1.XDataSource = 'xx';
SP1.YDataSource = 'yy';

%axis lims
axis equal
xlim([0 3]);
ylim([-1 2]);

%title
ETIME = 0;
timetitle = sprintf("%.3f seconds",ETIME);
title(timetitle);

%% Recording loop
i=1;
TMR=tic;
while 1
    buffer = deviceReader();
    
    CH.L = buffer(:,1);
    CH.R = buffer(:,2);
    
    AOA = [AOA; detect_az1(CH,CR,C,REC.d)];
    
%     if ~isplaying(droneObj)
%         play(droneObj);
%     end

    xx = [r.x r.x+d*cosd(90-AOA(i,1))];
    yy = [r.y r.y+d*sind(90-AOA(i,1))];
    
    %title
    ETIME = ETIME + BLK_t;
    timetitle = sprintf("%.3f seconds",ETIME);
    title(timetitle);
    
    refreshdata;
    drawnow;
    i=i+1;

end
clear TMR;

release(deviceReader)