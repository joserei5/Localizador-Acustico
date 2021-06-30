clear;clc;close all;
addpath ../../simLocUS
addpath ../../functions
addpath ../../soundfiles/generated
addpath ../../structures

%% Variable Definitions
% audio
WAV = 'mosquito_line_2.wav';
[s,f] = audioread(WAV);
Ns = length(s);

% room struct
room.xy = [5 7];
room.rc = [0 0 0 0];
room.T = 20;

% receiver struct
load('recstruct.mat');
rec.type = 0;
rec.loc = [1 3.5 1.7];
rec.mic.dmf = 1;
rec.mic.dist = 0.1*rec.mic.dmf*2;
rec.th = 0;
rec.phi = 0;

% plot points init
src.x=[];
src.y=[];

% trajectory
YAXIS=          1; % line across y-axis
XAXIS=          0; % line across x-axis
TRAJ_Y=         7; % trajectory y
TRAJ_X=         3; % trajectory x
TRAJ_VELOCITY=  1; % trajectory velocity

% block size
BLK_T=  100e-3; %100ms


%% Preparation of Values
% initialize Angle-of-arrival array (azimuth values)
% AOA = zeros(d.pts,1);

% split trajectory into blocks
BLK_S=      BLK_T*f;
BLK_INTERV= round(length(s)/BLK_S);

% trajectory(source) points
D_F=        TRAJ_VELOCITY/f;
SRC=        zeros(BLK_INTERV,BLK_INTERV);
if YAXIS
    SRC_Y=      0.1:D_F:(TRAJ_Y-0.1);
    SRC(:,2)=   SRC_Y(1:BLK_S:end-1);
    SRC(:,1)=   ones(BLK_INTERV,1)*TRAJ_X;
elseif XAXIS
    SRC_X=      0.1:D_F:(TRAJ_X-0.1);
    SRC(:,2)=   SRC_X(1:BLK_S:end-1);
    SRC(:,1)=   ones(BLK_INTERV,1)*TRAJ_Y;
end

% determine theoretical value of AOA
AOA_t = 180-atan2d(SRC(:,1)-rec.loc(1), SRC(:,2)-rec.loc(2));


%% Figure
% room
trjF = figure();
trjF.ToolBar = 'none';
trjF.MenuBar = 'none';
trjF.Resize = 'off';
rectangle('Position', [0 0 room.xy(1) room.xy(2)])
axis equal
xlim([0 room.xy(1)])
ylim([0 room.xy(2)])

% receiver
hold on;
rec_offset=rec.mic.pos(1,2);
plot(rec.loc(1), rec.loc(2), 'r+')
plot(rec.loc(1), rec.loc(2)+rec_offset, 'r>')
plot(rec.loc(1), rec.loc(2)-rec_offset, 'r>')
xlim([0 room.xy(1)])
ylim([0 room.xy(2)])
hold off;

% detection line
hold on;
DETx = [rec.loc(1) src.x];
DETy = [rec.loc(2) src.y];
DETp= plot(DETx, DETy, 'r');
DETp.XDataSource = 'DETx';
DETp.YDataSource = 'DETy';
hold off;

% theoretical azimuth line
hold on;
REALx = [rec.loc(1) src.x];
REALy = [rec.loc(2) src.y];
REALp= plot(REALx, REALy, 'k');
REALp.XDataSource = 'REALx';
REALp.YDataSource = 'REALy';
hold off;

%% Prepare blocks
tBLK = 100e-3;
szBLK = tBLK * f;
nBLK = ceil(Ns/szBLK);

% s_ = zeros(nBLK, szBLK);

% for i=1:nBLK
%     i1 = szBLK*(i-1)+1;
%    if i ~= nBLK
%        i2 = szBLK*i;
%        s_(i,:) = s(i1:i2,1);
%    else
%        i2 = i1+round((Ns/szBLK-(nBLK-1))*szBLK-1);
%        i3 = nBLK*szBLK - i2;
%        s_(i,:) = [s(i1:i2,1); zeros(i3,1)];
%    end
% end


%% Audio Device Writer -- init
fileReader = dsp.AudioFileReader(WAV);
fileReader.SamplesPerFrame = szBLK;

fileInfo = audioinfo(WAV)
deviceWriter = audioDeviceWriter('SampleRate',fileInfo.SampleRate);

bufferLatency = fileReader.SamplesPerFrame/deviceWriter.SampleRate

setup(deviceWriter,zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels))
info(deviceWriter)

CH.fs = fileInfo.SampleRate;
CR = round(rec.mic.dist + 1e-3*CH.fs);
C = 20.05*sqrt(273.15+room.T);
D_X = rec.mic.dist;

%% Audio Device Writer -- loop
i=1;
while ~isDone(fileReader)
    audioData = fileReader();
    
    CH.L = audioData(:,1);
    CH.R = audioData(:,2);
    
    AOA = detect_az1(CH,CR,C,D_X);
    
%     det_d = sqrt((rec.loc(1)-spk2.loc(1))^2+(rec.loc(2)-spk2.loc(2))^2);
    det_d = 3;
    DETx = [rec.loc(1)+det_d*cosd(AOA-90) rec.loc(1)];
    DETy = [rec.loc(2)+det_d*sind(AOA-90) rec.loc(2)];
    
%     REALx = [rec.loc(1)+det_d*cosd(AOA_t(i)-90) rec.loc(1)];
%     REALy = [rec.loc(2)+det_d*sind(AOA_t(i)-90) rec.loc(2)];
    i=i+1;
    
    refreshdata;
    drawnow;
    
    deviceWriter(audioData);
end

release(fileReader)
release(deviceWriter)