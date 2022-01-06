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
Nth = length(REC.th);

% [Audio settings]
AUDIO.name = 'mosquito2.wav';   % audio(source) file name
% AUDIO.name = 'drone1.mp3';      % audio(source) file name
BLK_t = 100*1e-3;               % block time size
N = 10/BLK_t;                   % number of blocks used
fl = 100;                       % lower frequency bound
fh = 18e3;                      % higher frequency bound

% [Noise settings]
SNR = -20:5:20;
Nsnr = length(SNR);
signalpower = 'measured';

% [Trajectory settings]
% -tm |----- 0 -----| tm 
tm = 1; % trajectory half-length 

% [Algorithm settings]
MED_t = 1e-3;                       % maximum estimated delay
C = 20.05*sqrt(273.15+ROOM.T);      % sound velocity
%--
AOA.theoretical.values = zeros(N,Nth*Nsnr);
%--
AOA.simulator.values = zeros(N,Nth*Nsnr);
AOA.simulator.error = zeros(N,Nth*Nsnr);
%--
AOA.algorithm.values = zeros(N,Nth*Nsnr);
AOA.algorithm.error = zeros(N,Nth*Nsnr);


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
errorf=figure;  % error figure (AOA errors)


%% Predictor Loop
% for azimuth=-45º and azimuth=+45º

idx=1; % global index
for k=1:Nsnr
    fprintf("SNR=%d dB\n", SNR(k));
    for i=1:Nth
        % [Trajectory]________________________________________
        TRAJ.y = linspace(-tm,tm,N)' + REC.xyz(2);
        TRAJ.x = ones(N,1)*(REC.xyz(1)+REC.d);
        TRAJ.z = ones(N,1)*REC.xyz(3);
        TRAJ.xyz = [TRAJ.x TRAJ.y TRAJ.z];

        % [Theoterical Results]______________________________
        SRC.traj.x = TRAJ.x;
        SRC.traj.y = TRAJ.y;
        ROOM.rec.x = REC.xyz(1);
        ROOM.rec.y = REC.xyz(2);
        ROOM.rec.dx = REC.DX;
        ROOM.rec.azimuth = REC.th(i);
        AOA_ = getTrajAOA(SRC,ROOM);
        AOA.theoretical.values(:,idx) = AOA_.theoretical;
        AOA.algorithm.values(:,idx) = AOA_.algorithm;
        AOA.algorithm.error(:,idx) = abs( AOA_.theoretical - AOA_.algorithm );

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
            CH.L = awgn(CH.L, SNR(k), signalpower);
            CH.R = awgn(CH.R, SNR(k), signalpower);
            % start timer
            TMR=tic;
            [AOA.simulator.values(j,idx),~] = detect_az3(CH, CR, C, REC.DX); %pred v3 (spline interp)
            % register lap (timer)
            ELAPSED=ELAPSED+toc(TMR);
            % calculate AOA error
            AOA.simulator.error(j,idx) = abs(AOA.simulator.values(j,idx) - AOA.theoretical.values(j,idx));
        end
        % mean of elapsed time for the N blocks
        ELAPSED = ELAPSED/N;

        % main loop statistics
        fprintf("========== RECEIVER THETA: %dº ==========\n", REC.th(i));
        fprintf("Predictor avg time: %.3f us\n", ELAPSED*1e6);
        fprintf("AOA: min=%.1fº max=%.1fº\n",min(AOA.simulator.values(:,idx)),max(AOA.simulator.values(:,idx)));
        
        % increment global index
        idx = idx + 1;
    end
    fprintf('\n');
end

%% Final Results
% reshape vectors
AOA.algorithm.values = reshape(AOA.algorithm.values,[N*2,Nsnr]);
AOA.algorithm.error = reshape(AOA.algorithm.error,[N*2,Nsnr]);
AOA.theoretical.values = reshape(AOA.theoretical.values,[N*2,Nsnr]);
AOA.simulator.values = reshape(AOA.simulator.values,[N*2,Nsnr]);
AOA.simulator.error = reshape(AOA.simulator.error,[N*2,Nsnr]);

% error figure
figure(errorf);
% determine Rows x Columns (RxC)
for i=2:6
    if mod(Nsnr,i) == 0
        RxC=i;
        break;
    elseif i==6
        RxC = 6;
    end
end
% plot
for i=1:Nsnr
    subplot(RxC,RxC,i)
        plot(AOA.theoretical.values(:,i), AOA.algorithm.error(:,i), 'g')
        hold on;
        plot(AOA.theoretical.values(:,i), AOA.simulator.error(:,i), 'k')
        hold off;
        xlim([0 180])
        title_string = sprintf('SNR = %d dB', SNR(i));
        title(title_string);
        xlabel('AOA Reference (º)');
        ylabel('AOA Error (º)');
end
% legend('Algorithm','Simulator')
% xlabel('AOA Reference (º)')
% ylabel('AOA Error (º)')