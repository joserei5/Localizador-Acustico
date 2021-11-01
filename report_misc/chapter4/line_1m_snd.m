clear all;clc;close force all;
addpath ../../functions
addpath ../../simLocUS

%% INIT
% Room
ROOM.xyz = [3 3 3]; % room coordinates
ROOM.T = 24; % room temperature
ROOM.H = 75; % room humidity
ROOM.P = 1.01; % room pressure
% Receiver
REC.xyz = [0.7 ROOM.xyz(2)/2 1.7]; % receiver coordinates
REC.d = 1; % receiver distance in relation to the source
REC.DX = 29.2*1e-2; % receiver: microphone interdistance
REC.th = -45; % receiver theta
REC.phi = 0; % receiver phi
% Audio settings
AUDIO.name = 'mosquito2.wav'; % audio(source) file name
AUDIO.fl = 100; % lower frequency bound
AUDIO.fh = 18e3; % higher frequency bound
% Trajectory
N = 9; % Number of points
tm = 1; % trajectory limit (+m and -m)
TRAJ.data.y = linspace(-tm,tm,N)' + REC.xyz(2);
TRAJ.data.x = ones(N,1)*(REC.xyz(1)+REC.d);
TRAJ.data.z = ones(N,1)*REC.xyz(3);
TRAJ.data.v = 0.2;

%% CREATE IDEAL ENVIRONMENT TRAJECTORY
ideal_wav(AUDIO, TRAJ, ROOM, REC);
%%
REC.th = 45; % receiver theta
ideal_wav(AUDIO, TRAJ, ROOM, REC);