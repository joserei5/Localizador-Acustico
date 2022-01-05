%% Include paths + Clean script
clear;clc;close all;
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../simLocUS
addpath ../../structures


%% Variables and Structures
% [Room]
ROOM.xyz = [3 3 3]; % coordinates
ROOM.T = 24;        % temperature
ROOM.H = 75;        % humidity

% [Receiver]
REC.xyz = [0.7 ROOM.xyz(2)/2 1.7];  % coordinates
REC.d = 1;                          % distance in relation to the source
REC.DX = 29.2*1e-2;                 % microphone interdistance
REC.th = [-45 45];                  % theta (a.k.a. azimuth)

% [Audio settings]
AUDIO.name = 'mosquito2.wav';   % audio(source) file name
BLK_t = 100*1e-3;               % block time size
N = 10/BLK_t;                   % number of blocks used
fl = 100;                       % lower frequency bound
fh = 18e3;                      % higher frequency bound

% [Trajectory settings]
% -tm |----- 0 -----| tm 
tm = 1; % trajectory half-length 

% [Algorithm settings]
MED_t = 1e-3;                       % maximum estimated delay
C = 20.05*sqrt(273.15+ROOM.T);      % sound velocity
AOA = zeros(N,length(REC.th));      % angle of arrival
AOA_error = AOA;                    % angle of arrival error
TRAJ.th = zeros(N,length(REC.th));  % trajectory theta/azimuth

AOA2.theoretical = zeros(N,length(REC.th));
AOA2.algorithm = zeros(N,length(REC.th));
AOA2.error = zeros(N,length(REC.th));


%% Process variables
% [Audio settings]
[AUDIO.y, CH.fs] = audioread(AUDIO.name); % read audio sample

% [Algorithm settings]
MED_N = MED_t * CH.fs;  % maximum estimated delay (samples)
CR = round(MED_N);      % correlation range
BLK_N = BLK_t * CH.fs;  % block size (samples)

% [sim_stereo.m]
% receiver structure
SS.REC.type = 0;                    % no surfaces around the microphone
SS.REC.struct = 'recstruct';        % default structure:
                                    % |_2Mic: 20cm diameter; dir: +yy
SS.REC.loc = REC.xyz;               % location in coordinates (x,y,z)
SS.REC.phi = 0;                     % 0º of elevation (normal)
SS.REC.mic.dmf = REC.DX/20e-2;      % distance multiplying factor:
                                    %  |_Usage: 2/diameter * new_diameter/2
                                    %           <=> new_diameter/diameter
                                    %  (1) [0 1 0]:
                                    %      [0 diameter/2 0] * 2/diameter
                                    %  (2) [0 new_diameter/2 0]:
                                    %      [0 1 0] * new_diameter/2
% room structure
SS.ROOM.size = ROOM.xyz;            % dimensions
SS.ROOM.coeff = [0 0 0];            % reflection coefficients:
                                    %  |_[walls ceiling floor]
SS.ROOM.MR = 0;                     % no reflections
SS.ROOM.R = Room();                 % create room object:
SS.ROOM.R.T = ROOM.T;               %  |_temperature
SS.ROOM.R.H = ROOM.H;               %  |_humidity
SS.ROOM.R.fs = CH.fs;               %  |_sampling frequency
SS.ROOM.R.fl = fl;                  %  |_lower freq. bound
SS.ROOM.R.fh = fh;                  %  |_upper freq. bound
% audio structure
SS.AUDIO.fs = CH.fs;                % sampling frequency
SS.AUDIO.fl = fl;                   % lower freq. bound
SS.AUDIO.fh = fh;                   % upper freq. bound


%% Initialize figure handles     
mainf=figure;   % main figure (studio + AOA curves)
errorf=figure;  % error figure (AOA errors)


%% Predictor Loop
% for azimuth=-45º and azimuth=+45º

for i=1:length(REC.th)
    % [Trajectory]________________________________________
    TRAJ.y = linspace(-tm,tm,N)' + REC.xyz(2);
    TRAJ.x = ones(N,1)*(REC.xyz(1)+REC.d);
    TRAJ.z = ones(N,1)*REC.xyz(3);
    TRAJ.xyz = [TRAJ.x TRAJ.y TRAJ.z];
%     TRAJ.th(:,i) = linspace(45,135,N) + REC.th(i);

    % [Theoterical Results]______________________________
    SRC.traj.x = TRAJ.x;
    SRC.traj.y = TRAJ.y;
    ROOM.rec.x = REC.xyz(1);
    ROOM.rec.y = REC.xyz(2);
    ROOM.rec.dx = REC.DX;
    ROOM.rec.azimuth = REC.th(i);
    AOA_ = getTrajAOA(SRC,ROOM);
    AOA2.theoretical(:,i) = AOA_.theoretical;
    AOA2.algorithm(:,i) = AOA_.algorithm;
    AOA2.error(:,i) = abs( AOA_.theoretical - AOA_.algorithm );
    
    % [Figures]___________________________________________
    figure(mainf);
    subplot(2,length(REC.th),i)
    % empty room
    r1=rectangle('Position',[0 0 ROOM.xyz(1:2)]);
    % shed
    r2=rectangle('Position',[0 REC.xyz(2)-1/2 1 1]);
    r2.FaceColor = [.9 .9 .9];
    % receiver mark
    hold on;
    plot(REC.xyz(1), REC.xyz(2), 'r*')
    % receiver auxiliary line
    rec_line = [0 -REC.DX/2 0; 0 REC.DX/2 0];
    rec_line = rotationMatrix(rec_line, -REC.th(i), 0);
    rec_line = rec_line + REC.xyz;
    plot(rec_line(:,1), rec_line(:,2), 'r');
    plot(rec_line(:,1), rec_line(:,2), 'r>');
    hold off;
    % trajectory line
    hold on;
    plot(TRAJ.xyz(:,1), TRAJ.xyz(:,2),'b')
    plot(TRAJ.xyz(1,1), TRAJ.xyz(1,2),'b*')
    hold off;
    % adjust axis
    axis equal
    xlim([0 ROOM.xyz(1)])
    ylim([0 ROOM.xyz(2)])

    % [Predictor]________________________________________
    % update theta for "sim_stereo.m"
    SS.REC.th = REC.th(i);

    % elapsed time (timer)
    ELAPSED=0;
    % predictor loop
    for j=1:N
        % indexing
        i1 = 1+(j-1)*BLK_N;
        i2 = j*BLK_N;
        % retrieve block samples
        SS.AUDIO.s = AUDIO.y(i1 : i2, 1);
        % update speaker location
        SPK.loc = TRAJ.xyz(j,:);
        % generate stereo samples (L+R)
        [CH.L, CH.R] = sim_stereo(SS.REC, SS.ROOM, SPK, SS.AUDIO, 1);
        % apply white gaussian noise
        % to each channel (separate noise for each channel)
%         snr=-20;
%         signalpower = 'measured';
%         CH.L = awgn(CH.L, snr, signalpower);
%         CH.R = awgn(CH.R, snr, signalpower);
        % start timer
        TMR=tic;
        [AOA(j,i),~] = detect_az3(CH, CR, C, REC.DX); %pred v3 (spline interp)
        % register lap (timer)
        ELAPSED=ELAPSED+toc(TMR);
        % calculate AOA error
%         AOA_error(j,i) = abs(AOA(j,i) - TRAJ.th(j,i));
        AOA_error(j,i) = abs(AOA(j,i) - AOA2.theoretical(j,i));
    end
    % mean of elapsed time for the N blocks
    ELAPSED = ELAPSED/N;
    
    % main loop statistics
    fprintf("========== RECEIVER THETA: %dº ==========\n", REC.th(i));
    fprintf("Predictor avg time: %.3f us\n", ELAPSED*1e6);
    fprintf("AOA: min=%.1fº max=%.1fº\n",min(AOA(:,i)),max(AOA(:,i)));
    
end


%% Experimental results
AOA3 = load('../../AudioInterface/experimental.mat');


%% Final Results
% Adjust labels
figure(mainf)
subplot(2,length(REC.th),1)
xlabel('m','Rotation',0)
ylabel('m','Rotation',0)
subplot(2,length(REC.th),2)
xlabel('m','Rotation',0)
ylabel('m','Rotation',0)

% reshape vectors
TRAJ.th = reshape(TRAJ.th,[],1);
AOA = reshape(AOA,[],1);
AOA_error = reshape(AOA_error,[],1);
AOA2.theoretical = reshape(AOA2.theoretical,[],1);
AOA2.algorithm = reshape(AOA2.algorithm,[],1);
AOA2.error = reshape(AOA2.error,[],1);

% Main figure
figure(mainf);
subplot(2,length(REC.th),length(REC.th)+1:2*length(REC.th));
plot(AOA2.theoretical, AOA2.theoretical, 'r:')
hold on;
plot(AOA2.theoretical, AOA2.algorithm, 'g')
plot(AOA2.theoretical, AOA,'k')
plot(AOA3.theoretical, AOA3.experimental.values,'b')
hold off;
xlabel('AOA Reference (º)')
ylabel('AOA Predictor (º)')
legend('Reference','Algorithm','Simulator','Experimental')

% Error figure
figure(errorf);
hold on;
plot(AOA2.theoretical, AOA2.error, 'g')
plot(AOA2.theoretical, AOA_error, 'k')
plot(AOA3.theoretical, AOA3.experimental.error,'b')
hold off;
legend('Algorithm','Simulator','Experimental')
xlabel('AOA Reference (º)')
ylabel('AOA Error (º)')