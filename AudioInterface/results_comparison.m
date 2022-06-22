clear;clc;close all;
load experimental;
addpath ../functions/

%% Variables
% [Room]
ROOM.xyz = [3 3 3]; % coordinates
ROOM.T = 24;        % temperature
ROOM.H = 75;        % humidity

% [Receiver]
REC.xyz = [0.7 ROOM.xyz(2)/2 1.7];  % coordinates
REC.d = 1;                          % distance in relation to the source
REC.DX = 29.2*1e-2;                 % microphone interdistance
REC.th = [-45 45];                  % theta (a.k.a. azimuth)
Nth = length(REC.th);

% [Temporal blocks]
BLK_t = 250*1e-3;               % block time size
N = 10/BLK_t;                   % number of blocks used

% [Trajectory settings]
% -tm |----- 0 -----| tm 
% trajectory half-length 
tm = 1; 

% [Algorithm settings]
CH.fs = 48e3;
C = 20.05*sqrt(273.15+ROOM.T);      % sound velocity

% [Structures]
AOA.theoretical.values = zeros(N*Nth,1);
AOA.algorithm.values = zeros(N*Nth,1);
AOA.algorithm.error = zeros(N*Nth,1);
AOA.simulator.values = zeros(N*Nth,1);
AOA.simulator.error = zeros(N*Nth,1);

%% Trajectory
% [Trajectory - w/o delay]
TRAJ.v = 0.2;
BLK_d = BLK_t * TRAJ.v;
TRAJ.y = transpose((REC.xyz(2)-tm)+BLK_d/2:BLK_d:(REC.xyz(2)+tm));
TRAJ.x = ones(N,1)*(REC.xyz(1)+REC.d);
TRAJ.z = ones(N,1)*REC.xyz(3);
TRAJ.xyz = [TRAJ.x TRAJ.y TRAJ.z];
% distances
TRAJ.d = sqrt((REC.xyz(1)-TRAJ.x).^2+(REC.xyz(2)-TRAJ.y).^2);

% [Delays]
% Sound wave delay
% (from RECEIVER to SOURCE)
TAU = TRAJ.d / C;

% [Trajectory - w/ delay]
% fix trajectory block center
BLK_d_adj = (BLK_t + TAU) * TRAJ.v;
TRAJ.y_adj = zeros(N,1);
TRAJ.y_adj(1,1) = (REC.xyz(2)-tm) + (BLK_d_adj(1)/2);
for i=2:N
    TRAJ.y_adj(i,1) = TRAJ.y_adj(i-1,1) + BLK_d_adj(i);
end
TRAJ.xyz_adj = [TRAJ.x TRAJ.y_adj TRAJ.z];

% [getTrajAOA.m]
% SRC (source) structure
SRC.traj.x = TRAJ.x;
SRC.traj.y = TRAJ.y;
% ROOM structure
ROOM.rec.x = REC.xyz(1);
ROOM.rec.y = REC.xyz(2);
ROOM.rec.dx = REC.DX;
% Delay estimation error (epsilon)
sample_error = 1/2;
epsilon = [-sample_error/CH.fs -sample_error/CH.fs];
% Process theoretical results 
for i=1:Nth
    % indexing
    i1 = N*(i-1) + 1;
    i2 = N*i;
    % update theta
    ROOM.rec.azimuth = REC.th(i);
    % Process trajectory w/ delay estimation
        % retrieve values from function
        AOA_ = getTrajAOAwError(SRC,ROOM,C,epsilon(i));
        % update structures
        AOA.theoretical.values(i1:i2,1) = AOA_.theoretical;
        AOA.algorithm.values(i1:i2,1) = AOA_.algorithm;
        % algorithm structure error
        AOA.algorithm.error(i1:i2,1) = abs( AOA_.theoretical - AOA_.algorithm );
end

%% Figures
errorf=figure;
plot(theoretical, experimental.error,'k:','LineWidth',1)
hold on;
plot(AOA.theoretical.values, AOA.algorithm.error, 'b')
hold off;

xlim([0 180])
xlabel('Aoa Reference (º)');
ylabel('AoA Error (º)');
legend('Expected results (w/ error)', 'Experiment results')