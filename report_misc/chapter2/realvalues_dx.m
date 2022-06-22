%% Initial cmds
clear;clc;close all;
addpath ../../simLocUS
addpath ../../functions

%% Data
delta_x = 29.2 * 1e-2;
dx_range = delta_x + linspace(1e-3,1e-2,50);
Ndx = length(dx_range);
room.temp = 20;
rcoord = [11 11 3]; % room coordinates
hcoord = [0.5 rcoord(2)/2 1.7]; % head coordinates
C = 20.05*sqrt(273.15+room.temp);

% trajectory
t_dist = 1;
t_points = 500;
t_offset = 0;
tcoord=speaker_hcircle2D(rcoord,hcoord,t_dist,t_points,t_offset);

% theoretical AOA and
% algorithm AOA
SRC.traj.x = tcoord(:,1);
SRC.traj.y = tcoord(:,2);
ROOM.rec.x = hcoord(1);
ROOM.rec.y = hcoord(2);
ROOM.rec.azimuth = 0;

% initialize figure
error_f = figure;

% initialize error arrays
max_error = zeros(Ndx,1);

for i=1:Ndx
    % update \Delta x
    ROOM.rec.dx = dx_range(i);
    % update the AOA values
    AOA = getTrajAOA(SRC,ROOM);
    % determine the AOA error
    AOA.error = abs(AOA.theoretical - AOA.algorithm);
    max_error(i,1) = max(AOA.error);
end

% transform the errors in percentages
percentage_error = abs(max_error(1)-max_error)*100;

% figure labels
plot(percentage_error)
xlabel('Theoretical AOA (º)')
ylabel('Predicted AOA error (º)')