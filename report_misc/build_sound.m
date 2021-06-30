clear;clc;close all;

addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generic
addpath ../structures
addpath ../structures/trajectories

%% read WAV audio
WAV=        'mosquito.wav';
[s, fs]=    audioread(WAV);
Ns=         length(s);
t=          (1:Ns)*1/fs;

% frequency limits
fl= 100;   % lower frequency bound
fh= 20e3;  % upper frequency bound (<=80 kHz due to @KemoL10_TF)

%% trajectory
load    custom150621_1110
TRAJ.X = TRAJ.X';
TRAJ.Y = TRAJ.Y';
NT = length(TRAJ.X);

% velocity
v=      0.25;

% trajectory jump
df=     v/fs;

% interpolate multiple curves
% F = scatteredInterpolant(TRAJ.X',TRAJ.Y', ones(length(TRAJ.X),1));
% d.x = F.Points(:,1);
% d.y = F.Points(:,2);
% clear F TRAJ

% rebuild trajectory
difX=diff(TRAJ.X);
difY=diff(TRAJ.Y);


d.x=TRAJ.X(1);
d.y=TRAJ.Y(1);

for i=2:NT
    % ALONG THE Y AXIS
   if       (difX(i)==0) && (difY(i)~=0)
       traj_done_y = TRAJ.Y(i-1):df:TRAJ.Y(i);
       d.y=[d.y traj_done_y];
       d.x=[d.x d.x(i-1)*ones(1,length(traj_done_y))];
    % ALONG THE X AXIS
   elseif   (difY(i)==0) && (difX(i)~=0)
       traj_done_x = TRAJ.X(i-1):df:TRAJ.X(i);
       d.x=[d.x traj_done_x];
       d.y=[d.y d.y(i-1)*ones(1,length(traj_done_x))];
    % ALONG BOTH AXIS
   else
       % STRONGER Y-AXIS or STRONGER X-AXIS
       if   (difY(i) > difX(i)) || (difY(i) < difX(i))
           Ns = length(TRAJ.X(i-1):df:TRAJ.X(i));
           vratio.x = (difX(i)/(difX(i)+difY(i)))*df;
           vratio.y = (difY(i)/(difX(i)+difY(i)))*df;
           traj_done_x = TRAJ.X(i-1):vratio.x:TRAJ.X(i);
           traj_done_y = TRAJ.Y(i-1):vratio.x:TRAJ.Y(i);
%            Ntd.x = length(traj_done_x);
%            Ntd.y = length(traj_done_y);
%            if Nt
            d.x=[d.x traj_done_x];
            d.y=[d.y traj_done_y];
           
       % EMPTY DISTANCES
       elseif (difY(i)==0) && (difX(i)==0)
           % do nothing
           d.x = [d.x TRAJ.X(i-1)];
           d.y = [d.y TRAJ.Y(i-1)];
           
       % DISTANCE OF X-AXIS = DISTANCE OF Y-AXIS
       else
           % SAME DISTANCES (HALF VELOCITY)
           traj_done_x = TRAJ.X(i-1):(df/2):TRAJ.X(i);
           traj_done_y = TRAJ.Y(i-1):(df/2):TRAJ.Y(i);
           d.x=[d.y traj_done_x];
           d.y=[d.y traj_done_y];
       end
   end
end


% N=      length(p_y);

% initialize stereo channels
% y=  zeros(2,N);