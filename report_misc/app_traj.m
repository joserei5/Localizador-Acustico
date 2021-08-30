clear all;clc;close all;
addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generated
addpath ../structures/proc_trajectories
addpath ../structures/trajectories


%% Parameters
% interp structure
interp = load('traj_interp040821_2159_49.mat');
% load('traj_interp010821_2242_22.mat');
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

% TRAJ structure
TRAJ = load('custom040821_2121_59.mat');
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

% sound sample
% snd = interp.snd.name;
% splitting to make sure we have the sample name only
snd = split(TRAJ.spatial.path,'/');
snd = snd{end};
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

% room properties
r.coord = interp.room.xy;
r.rc= zeros(1,3);
r.temp = 20;
r.hum = 55;
r.press = 1.01;
r.rorder = 0;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

% receiver properties
prec.loc = interp.rec.loc;
prec.dir= [0 0];
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

% detection parameters
det.x = [];
det.y = [];
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

% distance
src.x = interp.x;
src.y = interp.y;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬


%% Audio Setup
% load mosquito sample
[s, fs] = audioread(snd);
Ns = length(s);
t = (1:Ns)*1/fs;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬


%% Figure
% room
trjF = figure();
trjF.ToolBar = 'none';
trjF.MenuBar = 'none';
trjF.Resize = 'off';
trjF.Position = [trjF.Position(1:2) 1280   720];
rectangle('Position', [0 0 r.coord(1) r.coord(2)])

% shadow region
hold on
rectangle('Position', [0 0 prec.loc(1) r.coord(2)], 'FaceColor', [.9 .9 .9])
hold off

% trajectory
hold on;
plot(src.x(1), src.y(1), 'b+', 'LineWidth', 2)
plot(src.x, src.y, 'b')
plot(src.x(end), src.y(end), 'b*', 'LineWidth', 2)
hold off;

% receiver
hold on;
drec = 0.2;
recX = drec/2;
plot(prec.loc(1), prec.loc(2), 'r+')
plot(prec.loc(1), prec.loc(2)+recX, 'r>')
plot(prec.loc(1), prec.loc(2)-recX, 'r>')
hold off;

% detection line
hold on;
DETx = [prec.loc(1) det.x];
DETy = [prec.loc(2) det.y];
DETp= plot(DETx, DETy, 'r');
DETp.XDataSource = 'DETx';
DETp.YDataSource = 'DETy';
hold off;

% theoretical azimuth line
hold on;
% REALx = [prec.loc(1) src.x(1)];
% REALy = [prec.loc(2) src.y(2)];
REALx = src.x(1);
REALy = src.y(2);
REALp = plot(REALx, REALy, 'k*', 'LineWidth', 2);
REALp.XDataSource = 'REALx';
REALp.YDataSource = 'REALy';
hold off;

% adjust scale
axis equal
xlim([0 r.coord(1)])
ylim([0 r.coord(2)])
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬


%% Audio Device Writer -- init
% block settings
tBLK = 100e-3;
szBLK = tBLK * fs;

% determine theoretical value of AOA
AOA_t = 180-atan2d(src.x-prec.loc(1), src.y-prec.loc(2));

% figure;
% plot(AOA_t)

% determine theoretical distances
% d_th = sqrt((prec.loc(1)-src.x).^2 + (prec.loc(2)-src.y).^2);
d = sqrt(r.coord(1)^2+r.coord(2).^2);

% figure;
% plot(d_th)

fileReader = dsp.AudioFileReader(snd);
fileReader.SamplesPerFrame = szBLK;

fileInfo = audioinfo(snd)
deviceWriter = audioDeviceWriter('SampleRate',fileInfo.SampleRate);

bufferLatency = fileReader.SamplesPerFrame/deviceWriter.SampleRate

setup(deviceWriter,zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels))
info(deviceWriter);

CH.fs = fileInfo.SampleRate;
CR = round(drec + 1e-3*CH.fs);
C = 20.05*sqrt(273.15+r.temp);
D_X = drec;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬


%% Audio Device Writer -- loop

i=1;
while ~isDone(fileReader)
    audioData = fileReader();
    
    CH.L = audioData(:,1);
    CH.R = audioData(:,2);
    
    AOA = detect_az1(CH,CR,C,D_X);
    
%     DETx = [prec.loc(1)+d_th(i)*cosd(AOA-90) prec.loc(1)];
%     DETy = [prec.loc(2)+d_th(i)*sind(AOA-90) prec.loc(2)];

    DETx = [prec.loc(1)+d*cosd(AOA-90) prec.loc(1)];
    DETy = [prec.loc(2)+d*sind(AOA-90) prec.loc(2)];
    
%     REALx = [prec.loc(1)+d_th(i)*cosd(AOA_t(i)-90) prec.loc(1)];
%     REALy = [prec.loc(2)+d_th(i)*sind(AOA_t(i)-90) prec.loc(2)];
    REALx = src.x(i);
    REALy = src.y(i);
    i=i+(4800-1);
    
%     [i cosd(AOA_t(i)) sind(AOA_t(i))]
    
    refreshdata;
    drawnow;
    
    deviceWriter(audioData);
end

release(fileReader)
release(deviceWriter)