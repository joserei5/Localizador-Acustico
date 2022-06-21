%% Include paths + Clean script
clear;clc;close force all;
addpath ../../functions
addpath ../../soundfiles/generated
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
AUDIO.name = {  'neg45_line6.wav'...  % audio(source) file name
                'pos45_line6.wav'   };
    
BLK_t = 250*1e-3;               % block time size
N = 10/BLK_t;                   % number of blocks used
fl = 100;                       % lower frequency bound
fh = 18e3;                      % higher frequency bound

% [Noise settings]
SNR = [-100 -50 -30 -19 -16 -13 -10 0 20]; % signal to noise ratio (dB)
Nsnr = length(SNR);                        % length of SNR array
signalpower = 'measured';                  % string shortcut

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
REP = 100000;


%% Process variables
% [Algorithm settings]
CH.fs = 48e3;           % sampling frequency is 48kHz (default)
MED_N = MED_t * CH.fs;  % maximum estimated delay (samples)
CR = round(MED_N);      % correlation range
BLK_N = BLK_t * CH.fs;  % block size (samples)

% [Trajectory - w/o delay]
TRAJ.v = 0.2;                                       % traj. velocity
BLK_d = BLK_t * TRAJ.v;                             % size (m) of a block
                                                    % in a line traj.
                                                    %
                                                    % traj. y axis:
TRAJ.y = transpose( (REC.xyz(2)-tm)+BLK_d/2 :...    % Start: middle of block;
                    BLK_d                   :...    % Increment: size of block;
                    (REC.xyz(2)+tm));               % End: middle of last block.
                                                    %
TRAJ.x = ones(N,1)*(REC.xyz(1)+REC.d);              % traj. x axis
TRAJ.z = ones(N,1)*REC.xyz(3);                      % traj. z axis
TRAJ.xyz = [TRAJ.x TRAJ.y TRAJ.z];                  % all axis traj. data

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

% initialize clock
CLK = tic;

% Experience repetitions loop
for NE=1:REP
    
    % global index
    % SAME AS
    % column counter;
    % TIP: it resets when another experience is made, i.e.,
    % when ALL the SNR values have been processed.
    idx=1;

    % SNR vector sweep loop
    for k=1:Nsnr
        
        % vector index
        % SAME AS
        % row counter;
        % TIP: it resets when the column changes, i.e.,
        % when the SNR changes.
        vidx=1;
        
        % Theta vector sweep loop
        for i=1:Nth
            
            % Load respective audio sample
            [AUDIO.y, CH.fs] = audioread(AUDIO.name{i}); % read audio sample
            
            % Predictor loop
            for j=1:N
                % indexing
                i1 = 1+(j-1)*BLK_N;
                i2 = j*BLK_N;
                
                % update speaker location
                SPK.loc = TRAJ.xyz(j,:);
                
                % retrieve stereo samples (L+R)
                CH.L = AUDIO.y(i1:i2 , 1);
                CH.R = AUDIO.y(i1:i2 , 2);
                
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
                % OR
                % change to the next row
                vidx = vidx + 1;
            end
        end
        
        % increment global index
        % OR
        % change to the next column
        idx = idx + 1;
    end
    
    fprintf("%d\n",NE);

end

% average the values
AOA.simulator.values = AOA.simulator.values/REP;
AOA.simulator.error = AOA.simulator.error/REP;

% register time
TMR = toc(CLK);
clear CLK;
% print time
TMR_.h = floor(TMR/3600);
TMR_.m = floor(mod(TMR,3600)/60);
TMR_.s = floor(mod(mod(TMR,3600),60));
fprintf("%2dh%2dm%2ds\n",TMR_.h, TMR_.m, TMR_.s);

% close waitbar
% close(wb);


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
for i=4:Nsnr
    subplot(2,RxC,i-3)
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