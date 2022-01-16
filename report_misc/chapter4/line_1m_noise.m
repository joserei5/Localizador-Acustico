%% Include paths + Clean script
clear;clc;close force all;
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
% AUDIO.name = 'mosquito2.wav';   % audio(source) file name
AUDIO.name = 'drone1.mp3';      % audio(source) file name
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
% trajectory half-length 
tm = 1; 

% [Algorithm settings]
MED_t = 1e-3;                       % maximum estimated delay
C = 20.05*sqrt(273.15+ROOM.T);      % sound velocity

% [Structures]
AOA.theoretical.values = zeros(N*Nth,1);
AOA.algorithm.values = zeros(N*Nth,1);
AOA.algorithm.error = zeros(N*Nth,1);
AOA.simulator.values = zeros(N*Nth,Nsnr);
AOA.simulator.error = zeros(N*Nth,Nsnr);


% [Experience variables]
% Number of experiences/repetitions
REP = 1000;


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

% [Room structure]
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

% [Audio structure]
SS.AUDIO.fs = CH.fs;                % sampling frequency
SS.AUDIO.fl = fl;                   % lower freq. bound
SS.AUDIO.fh = fh;                   % upper freq. bound

% [Trajectory]
TRAJ.y = linspace(-tm,tm,N)' + REC.xyz(2);
TRAJ.x = ones(N,1)*(REC.xyz(1)+REC.d);
TRAJ.z = ones(N,1)*REC.xyz(3);
TRAJ.xyz = [TRAJ.x TRAJ.y TRAJ.z];

% [getTrajAOA.m]
% SRC (source) structure
SRC.traj.x = TRAJ.x;
SRC.traj.y = TRAJ.y;
% ROOM structure
ROOM.rec.x = REC.xyz(1);
ROOM.rec.y = REC.xyz(2);
ROOM.rec.dx = REC.DX;
for i=1:Nth
    % indexing
    i1 = N*(i-1) + 1;
    i2 = N*i;
    % update theta
    ROOM.rec.azimuth = REC.th(i);
    % retrieve values from function
    AOA_ = getTrajAOA(SRC,ROOM);
    % update structures
    AOA.theoretical.values(i1:i2,1) = AOA_.theoretical;
    AOA.algorithm.values(i1:i2,1) = AOA_.algorithm;
    % algorithm structure error
    AOA.algorithm.error(i1:i2,1) = abs( AOA_.theoretical - AOA_.algorithm );
end



%% Predictor Loop
% for azimuth=-45º and azimuth=+45º

% Initialize waitbar
wb = waitbar(0,['0/' num2str(REP)]);

% initialize clock
CLK = tic;

% Experience repetitions loop
for NE=1:REP
    
    % global index
    idx=1;

    % SNR vector sweep loop
    for k=1:Nsnr
        
        % vector index
        vidx=1;
        
        % Theta vector sweep loop
        for i=1:Nth
            % Predictor loop
            for j=1:N
                % indexing
                i1 = 1+(j-1)*BLK_N;
                i2 = j*BLK_N;
                
                % retrieve block samples
                SS.AUDIO.s = AUDIO.y(i1 : i2, 1);
                
                % update speaker location
                SPK.loc = TRAJ.xyz(j,:);
                
                % update SS structure's theta
                SS.REC.th = REC.th(i);
                % generate stereo samples (L+R)
                [CH.L, CH.R] = sim_stereo(SS.REC, SS.ROOM, SPK, SS.AUDIO, 1);
                
                % apply white gaussian noise
                % to each channel (separate noise for each channel)
                CH.L = awgn(CH.L, SNR(k), signalpower);
                CH.R = awgn(CH.R, SNR(k), signalpower);
                
                % predictor v3 (spline interpolation)
                [EXP_AOA,~] = detect_az3(CH, CR, C, REC.DX);
                EXP_ERROR = abs(EXP_AOA - AOA.theoretical.values(vidx,1));
                
                % sum current values to previous values
                AOA.simulator.values(vidx,idx) = AOA.simulator.values(vidx,idx) + EXP_AOA;
                AOA.simulator.error(vidx,idx) = AOA.simulator.error(vidx,idx)+ EXP_ERROR;
                
                % increment vector index
                vidx = vidx + 1;
            end
        end
        % increment global index
        idx = idx + 1;
    end

    % update waitbar
    waitbar(NE/REP,wb,[num2str(NE) '/' num2str(REP)]);
end

% average the values
AOA.simulator.values = AOA.simulator.values/REP;
AOA.simulator.error = AOA.simulator.error/REP;

% register time
TMR = toc(CLK);
clear CLK;

% close waitbar
close(wb);


%% Final Results
% error figure
errorf=figure;
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
        plot(AOA.theoretical.values(:,1), AOA.algorithm.error(:,1), 'g')
        hold on;
        plot(AOA.theoretical.values(:,1), AOA.simulator.error(:,i), 'k')
        hold off;
        xlim([0 180])
        title_string = sprintf('SNR = %d dB', SNR(i));
        title(title_string);
        xlabel('AOA Reference (º)');
        ylabel('AOA Error (º)');
end