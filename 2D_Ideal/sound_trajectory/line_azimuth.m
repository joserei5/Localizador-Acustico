clear;clc;close all;
addpath ../../simLocUS
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../structures

%% Parameters
%▬SOUND FILE▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
IN_WAV= 'mosquito_line_2.wav';
%▬MAXIMUM ESTIMATED DELAY▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
MED=    1e-3;
%▬ROOM TEMPERATURE▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
R.T=    20;
%▬TRAJECTORY CHARACTERISTICS▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
YAXIS=          1; % line across y-axis
XAXIS=          0; % line across x-axis
TRAJ_Y=         7; % trajectory y
TRAJ_X=         3; % trajectory x
TRAJ_VELOCITY=  1; % trajectory velocity
%▬MICROPHONE CHARACTERISTICS▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
M_X=    1;
M_Y=    3.5;
D_X=    0.1;
%▬BLOCK SIZE▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
BLK_T=  100e-3; %100ms

%% Processing Parameters
% read audio
[s, fs]=    audioread(IN_WAV);
% sound velocity
C=          20.05*sqrt(273.15 + R.T);
% associate MED with CORRELATION RANGE
CORR_RANGE= floor((0.2/C)*fs + MED*fs);
% split audio into blocks
BLK_S=      BLK_T*fs;
BLK_INTERV= floor(length(s)/BLK_S);
% intialize azimuth array
DET_AZ=     zeros(BLK_INTERV,1);
% trajectory(source) points
D_F=        TRAJ_VELOCITY/fs;
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

%% Azimuth detection
for i=1:BLK_INTERV
    % indexes
    START_i=    (i-1)*BLK_S+1;
    END_i=      BLK_S*i;
    % channels
    L_CH=       s(START_i:END_i,1);
    R_CH=       s(START_i:END_i,2);
    % detect azimuth (algorithm)
    CH_CORR=            xcorr(L_CH,R_CH,CORR_RANGE);
    [CORR_V, CORR_i]=   max(abs(CH_CORR));
    DELAY_i=            (CORR_i-1) - CORR_RANGE;
    DELAY_t=            DELAY_i/fs;
    arg_=               C*DELAY_t/(2*D_X);
    arg_(arg_>1)=       1;
    arg_(arg_<-1)=      -1;
    DET_AZ(i)=          acosd(arg_);
end

%% Visualize detected values
figure;

SRC_X=  SRC(1,1);
SRC_Y=  SRC(1,2);
splt=   plot(SRC_X, SRC_Y, 'k+');
splt.XDataSource = 'SRC_X';
splt.YDataSource = 'SRC_Y';

pbaspect([1 1 1]) % axis square (quero é o axis equal)
xlim([0 5]);
ylim([0 7]);
hold on

L_X=    [M_X SRC_X];
L_Y=    [M_Y SRC_Y];
splt2=  plot(L_X, L_Y, 'k-+');
splt2.XDataSource = 'L_X';
splt2.YDataSource = 'L_Y';

AZ_X=   M_X+4*cosd(DET_AZ(1)-90);
AZ_Y=   M_Y+4*sind(DET_AZ(1)-90);
P_X=    [M_X AZ_X];
P_Y=    [M_Y AZ_Y];

aplt= plot(P_X, P_Y, 'r');
aplt.XDataSource = 'P_X';
aplt.YDataSource = 'P_Y';

D_DIFF=     zeros(BLK_INTERV,1);
AZ_DIFF=    zeros(BLK_INTERV,1);
for i=1:BLK_INTERV
    SRC_X=  SRC(i,1);
    SRC_Y=  SRC(i,2);
    d=      sqrt((M_X-SRC_X)^2+(M_Y-SRC_Y)^2);
    AZ_X=   M_X+d*cosd(DET_AZ(i)-90);
    AZ_Y=   M_Y+d*sind(DET_AZ(i)-90);
    P_X=    [M_X AZ_X];
    P_Y=    [M_Y AZ_Y];
    L_X=    [M_X SRC_X];
    L_Y=    [M_Y SRC_Y];
    refreshdata
    drawnow
    D_DIFF(i,1)=    sqrt((AZ_X-SRC_X)^2+(AZ_Y-SRC_Y)^2);
    AZ_DIFF(i,1)=   atand((SRC_Y-M_Y)/(SRC_X-M_X))...
                    - atand((AZ_Y-M_Y)/(AZ_X-M_X));
    
    pause(0.1)
end

figure;
subplot(211)
    plot(SRC(:,2),D_DIFF,'b')
    hold on
    plot(M_Y,0,'r+', 'LineWidth', 3)
    hold off
    xlabel('Trajectory (m)')
    ylabel('Distance from real source (m)')
subplot(212)
    plot(SRC(:,2),AZ_DIFF,'b')
    xlabel('Trajectory (m)')
    ylabel('[REAL] Azimuth error (º)')