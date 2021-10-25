clear all;clc;close all;
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../simLocUS
addpath ../../structures

%% Variables (Tune-in)
ROOM.xyz = [3 3 3]; % room coordinates
ROOM.T = 24; % room temperature
ROOM.H = 75; % room humidity
REC.xyz = [0.7 1.5 1.7]; % receiver coordinates
REC.d = 1; % receiver distance in relation to the source
REC.DX = 29.2*1e-2; % receiver: microphone interdistance

AUDIO.name = 'mosquito2.wav'; % audio(source) file name
N = 100; % number of blocks used
BLK_t = 100*1e-3; % block time size
fl = 100; % lower frequency bound
fh = 18e3; % higher frequency bound

MED_t = 1e-3; % maximum estimated delay
C = 20.05*sqrt(273.15+ROOM.T); % sound velocity
AOA = zeros(N,1); % angle of arrival

%% Create shed and room
fig=figure;
% empty room
r1=rectangle('Position',[0 0 ROOM.xyz(1:2)]);
% shed
r2=rectangle('Position',[0 1 1 1]);
r2.FaceColor = [.9 .9 .9];
% receiver mark
hold on;
plot(REC.xyz(1), REC.xyz(2), 'r+')
hold off;
% adjust axis
axis equal
xlim([0 ROOM.xyz(1)])
ylim([0 ROOM.xyz(2)])

%% Create Semi-circle
[TRAJ.xyz, TRAJ.th] = speaker_hcircle2D(ROOM.xyz, REC.xyz, REC.d, N, 0);
% [TRAJ.xyz, TRAJ.th] = speaker_circle2D(ROOM.xyz, REC.xyz+[.5 0 0], REC.d, N, 0);
% draw semi-circle
hold on;
plot(TRAJ.xyz(:,1), TRAJ.xyz(:,2),'b')
hold off;

%% Algorithm
% read audio
[AUDIO.y, CH.fs] = audioread(AUDIO.name);

% correlation range estimation
MED_N = MED_t * CH.fs;
CR = round(MED_N);
% get block sample size
BLK_N = BLK_t * CH.fs; % get audio samples <=> 100ms

% == SIM STEREO CONFIG ==
SS.REC.type = 0; % no surfaces around the microphone
SS.REC.struct = 'recstruct';
SS.REC.loc = REC.xyz;
SS.REC.th = 0;
SS.REC.phi = 0;
SS.REC.mic.dmf = 1;

SS.ROOM.size = ROOM.xyz;
SS.ROOM.coeff = [0 0 0];
SS.ROOM.MR = 0;
SS.ROOM.R = Room();
SS.ROOM.R.T = ROOM.T;
SS.ROOM.R.H = ROOM.H;
SS.ROOM.R.fs = CH.fs;
SS.ROOM.R.fl = fl;
SS.ROOM.R.fh = fh;

SS.AUDIO.fs = CH.fs;
SS.AUDIO.fl = fl;
SS.AUDIO.fh = fh;
% == EndOfConfig ==

% predictor loop
ELAPSED=0;
for i=1:N
    SS.AUDIO.s = AUDIO.y(1+(i-1)*BLK_N : i*BLK_N, 1);
    SPK.loc = TRAJ.xyz(i,:);
    [CH.L, CH.R] = sim_stereo(SS.REC, SS.ROOM, SPK, SS.AUDIO, 1);
    TMR=tic;
    [AOA(i,1),~] = detect_az3(CH, CR, C, REC.DX); %pred v3 (spline interp)
    ELAPSED=ELAPSED+toc(TMR);
%     [AOA(i,1),~] = detect_az2(CH, CR, C, REC.DX); %pred v2 (3point interp)
%     AOA(i,1) = detect_az1(CH, CR, C, REC.DX); %pred v1 (non-interp)
end
ELAPSED = ELAPSED/N;
fprintf("Mean time for predictor algorithm execution: %.3f us\n", ELAPSED*1e6);
fprintf("AOA: min=%.1fº max=%.1fº\n",min(AOA),max(AOA))

AOA_error = abs(TRAJ.th - AOA');

%% Final figures
figure;
plot(TRAJ.th, AOA)
hold on;
plot(TRAJ.th, TRAJ.th, 'r:')
hold off;

figure;
plot(TRAJ.th, AOA_error, 'b')