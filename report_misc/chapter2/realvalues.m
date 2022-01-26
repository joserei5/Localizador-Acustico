%% Initial cmds
clear;clc;close all;
addpath ../../simLocUS
addpath ../../functions

%% Data
delta_x = 29.2 * 1e-2;
room.temp = 20;
rcoord = [11 11 3]; % room coordinates
hcoord = [0.5 rcoord(2)/2 1.7]; % head coordinates
C = 20.05*sqrt(273.15+room.temp);

% trajectory
t_dist = 0.2;
t_points = 500;
t_offset = 0;
tcoord=speaker_hcircle2D(rcoord,hcoord,t_dist,t_points,t_offset);

% theoretical AOA and
% algorithm AOA
SRC.traj.x = tcoord(:,1);
SRC.traj.y = tcoord(:,2);
ROOM.rec.x = hcoord(1);
ROOM.rec.y = hcoord(2);
ROOM.rec.dx = delta_x;
ROOM.rec.azimuth = 0;
AOA = getTrajAOA(SRC,ROOM);

% AOA error
AOA.error = abs(AOA.theoretical - AOA.algorithm);

%% Plot
% room
roomf=figure;
rectangle('Position',[0 0 rcoord(1) rcoord(2)])
% trajectory
hold on;
plot(tcoord(:,1), tcoord(:,2),'k');
plot(tcoord(1,1), tcoord(1,2),'ko');
hold off;
% receiver
hold on;
plot(hcoord(1), hcoord(2), 'b+', 'LineWidth', 2);
hold off;

% results
figure;
plot(AOA.theoretical, AOA.theoretical, 'k');
hold on
plot(AOA.theoretical, AOA.algorithm, 'r');
hold off
legend('Theoretical AOA', 'Algorithm AOA')
xlabel('Theoretical AOA (º)')
ylabel('Predicted AOA (º)')

% errors
figure;
plot(AOA.theoretical, AOA.error, 'b')
xlabel('Theoretical AOA (º)')
ylabel('Predicted AOA error (º)')