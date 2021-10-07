clear force all;
clc;
% close force all;

%% Audio device info
qc.fs = 48e3;   % sampling frequency
qc.b  = 24;      % audio bits

%% Algorithm Variables
MED_t = 2e-3;                   % maximum estimated delay = 1 ms
MED_N = MED_t * qc.fs;          % maximum estimated delay (in samples)

ROOM.T = 23;                    % room temperature (ºC)
C = 20.05*sqrt(273.15 + ROOM.T);% sound velocity (m/s)
REC.d = 29.5e-2;                % receiver distance = 29.5 cm
% REC.dN = (REC.d/C)*qc.fs;       % receiver distance (in samples)

CH.fs = qc.fs;

% CR = round(REC.dN + MED_N);     % correlation range=receiver distance + MED (in samples)
CR = round(MED_N);     % correlation range=receiver distance + MED (in samples)

%% AudioWriter
snd = '../soundfiles/capture/0deg/recording290921_000412.wav'; % recording directory
BLK_t = 250e-3;             % Block size = 100 ms
BLK_N = BLK_t * qc.fs;      % Block size (in samples)

fileReader = dsp.AudioFileReader(snd);
fileReader.SamplesPerFrame = BLK_N;

fileInfo = audioinfo(snd)
deviceWriter = audioDeviceWriter('SampleRate',fileInfo.SampleRate);
bufferLatency = fileReader.SamplesPerFrame/deviceWriter.SampleRate

setup(deviceWriter,zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels))
info(deviceWriter)

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

%% Detection Loop
addpath ../functions
AOA = zeros(fileInfo.TotalSamples/BLK_N,1);

i=1;
while ~isDone(fileReader)
    audioData = fileReader();
    
    CH.L = audioData(:,1);
    CH.R = audioData(:,2);
    
    AOA(i,1) = detect_az1(CH,CR,C,REC.d);
    
    xx = [r.x r.x+d*cosd(90-AOA(i,1))];
    yy = [r.y r.y+d*sind(90-AOA(i,1))];
    
    %title
    ETIME = ETIME + BLK_t;
    timetitle = sprintf("%.3f seconds",ETIME);
    title(timetitle);
    
    refreshdata;
    drawnow;
    i=i+1;
    
    deviceWriter(audioData);
end

hold off;
release(fileReader)
release(deviceWriter)

%% Final AOA plot
figure;
t=BLK_t:BLK_t:10;
plot(t,AOA);
xlabel('seconds')
ylabel('AOA (degrees)')
title('AOA(t)')