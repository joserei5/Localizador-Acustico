%% Initial cmds
clear;clc;close all;
addpath ../../functions

%% TUNE-IN - 1 - POINTS
Nx = 100;       % no. points for Microphone Spacing (d_x)
Nd =  80;       % no. points for distance (d)

%% TUNE-IN - 2 - AMBIENT
room.temp = 20; % room temperature (ºC)

%% TUNE-IN - 3 - LIMITS
theta.min = 0;      % minimum AOA value
theta.max = 180;    % maximum AOA value
theta.steps = 15;   % AOA increment/steps

%% PROCESS VARIABLES - 1 - STATIC
C = 20.05*sqrt(273.15+room.temp);   % sound velocity

AOA = theta.min:theta.steps:theta.max;   % Nsnr range

d_x = linspace(0,2,Nx);                     % microphone spacing - range
d = linspace(0,7,Nd);                       % distance to the source - range

recL = [-d_x/2; zeros(1,Nx)];   % left microphone @ receiver
recR = [+d_x/2; zeros(1,Nx)];   % right microphone @ receiver

realAz = zeros(Nd,Nx);      % ref for \epsilon=0
azimuth = zeros(Nd,Nx);     % actual values
azimuth_re = ones(Nd,Nx);   % azimuth real error

%% [CICLE THROUGH ALL THETAS]
Nloop = length(AOA);
min_d1 = zeros(Nloop,1);
min_d5 = zeros(Nloop,1);

epsilon = 21 * 1e-6;

for i=1:Nloop
    %% PROCESS VARIABLES - 2 - UPDATEABLE
    actual_theta = AOA(i);                                                  % current theta or AOA
    source = [d'*cosd(actual_theta) d'*sind(actual_theta)];                 % source/mosquito location
    d1 = sqrt((recL(1, :)-source(:, 1)).^2 + (recL(2, :)-source(:, 2)).^2); % distance d1
    d2 = sqrt((recR(1, :)-source(:, 1)).^2 + (recR(2, :)-source(:, 2)).^2); % distance d2
    d_d = d1 - d2;                                                          % distance difference (m)
    d_t = d_d/C;                                                            % TDoA / time delay (s)

    %% AOA ALGORITHM / PREDICTOR
    arg = C*(d_t+epsilon)./d_x;
    arg(arg<-1) = -1;
    arg(arg>1) = 1;
    azimuth = acosd(arg);
    azimuth_re = abs(ones(Nd,Nx)*actual_theta - azimuth);

    %% RETRIEVE CONTOUR DATA
    % 1 degree error
    P1 = contourc(d_x, d, azimuth_re, [1 1]);
    % 5 degrees error
    P5 = contourc(d_x, d, azimuth_re, [5 5]);
    
    %% STORE THE MINIMAL DISTANCE REQUIRED
    if isempty(P1)
        min_d1(i,1) = 0;
    else
        min_d1(i,1) = max(P1(2,round(end/2):end))/2;
    end
    
    if isempty(P5)
        min_d5(i,1) = 0;
    else
        min_d5(i,1) = max(P5(2,round(end/2):end))/2;
    end
    
end

%% PLOT RESULTS
error_f = figure;
plot(AOA,min_d1)
hold on;
plot(AOA,min_d5)
hold off;
legend('1º error', '5º error')