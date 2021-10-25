clear;clc;close all;

addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generic
addpath ../structures
addpath ../structures/trajectories

%% read WAV audio
tic
WAV=        'mosquito.wav';
[s, fs]=    audioread(WAV); 
Ns=         length(s);
t=          (1:Ns)*1/fs;

% frequency limits
fl= 100;   % lower frequency bound
fh= 20e3;  % upper frequency bound (<=80 kHz due to @KemoL10_TF)

%% trajectory
TRAJ=load('VariableVelocity.mat');
THR = length(TRAJ.data.x);

waypoints = [TRAJ.data.x(1:THR), TRAJ.data.y(1:THR), TRAJ.data.z(1:THR)];

TRAJ.DX = diff(TRAJ.data.x(1:THR));
TRAJ.DX = [0; TRAJ.DX];
TRAJ.DY = diff(TRAJ.data.y(1:THR));
TRAJ.DY = [0; TRAJ.DY];

TRAJ.D = sqrt(TRAJ.DX.^2 + TRAJ.DY.^2);
% fill empty distances
% TRAJ.D(TRAJ.D == 0) = 1e-9;

TRAJ.T = TRAJ.D ./ TRAJ.data.v(1:THR);
TRAJ.T = cumsum(TRAJ.T);

s_rate = 48e3;
s_frames = s_rate;
trajectory = waypointTrajectory(waypoints,'TimeOfArrival',TRAJ.T,'SampleRate', s_rate, 'SamplesPerFrame', s_frames);

tx = [];
ty = [];
last_dump = round(round(TRAJ.T(end))*(s_rate/s_frames));
prev_disp=-1;
for i=1:last_dump
   dump = trajectory();
   tx = [tx; dump(:,1)];
   ty = [ty; dump(:,2)];
   
   idx=round(i/last_dump*100);
   if prev_disp ~= idx
       prev_disp = idx;
       fprintf("\r%d\n\b",idx)
   end
end


figure; plot(TRAJ.X, TRAJ.Y, 'b')
hold on;
plot(tx, ty, 'r');
hold off;
toc

% N=      length(p_y);

% initialize stereo channels
% y=  zeros(2,N);