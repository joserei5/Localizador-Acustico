clear all;clc;close all;
addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generic
addpath ../structures
addpath ../structures/trajectories

%%
VEL = 0.25;
TOA = (0:2:4) / VEL;
waypoints = [   0.1118 0.1880;  ...
                0.4827 0.8324;  ...
                0.9366 0.3688   ];
waypoints = [waypoints zeros(3,1)];

plot(waypoints(:,1), waypoints(:,2));

% SRATE = ceil(1000/(TOA(end)/3));
SRATE = 48e3;
FRAMES = SRATE;
trajectory = waypointTrajectory(waypoints,'TimeOfArrival',TOA, 'SampleRate', SRATE, 'SamplesPerFrame', FRAMES);

tx = [];
ty = [];
last_dump = round(round(TOA(end))*(SRATE/FRAMES));

% 349.73s /FULL
% 321.49  /HALF
tic
for i=1:last_dump
   dump = trajectory();
   tx = [tx; dump(:,1)];
   ty = [ty; dump(:,2)];
end
toc

hold on
plot(tx,ty)
hold off