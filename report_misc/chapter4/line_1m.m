clear all;clc;close all;
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../simLocUS
addpath ../../structures

%% Variables
% Room
ROOM.xyz = [3 3 3]; % room coordinates
ROOM.T = 24; % room temperature
ROOM.H = 75; % room humidity
% Receiver
REC.xyz = [0.7 ROOM.xyz(2)/2 1.7]; % receiver coordinates
REC.d = 1; % receiver distance in relation to the source
REC.DX = 29.2*1e-2; % receiver: microphone interdistance
REC.th = [-45 45]; % receiver theta
% Audio settings
AUDIO.name = 'mosquito2.wav'; % audio(source) file name
BLK_t = 100*1e-3; % block time size
N = 10/BLK_t; % number of blocks used
fl = 100; % lower frequency bound
fh = 18e3; % higher frequency bound
% Trajectory settings
tm = 1; % trajectory distance from center (90degrees)
% Algorithm settings
MED_t = 1e-3; % maximum estimated delay
C = 20.05*sqrt(273.15+ROOM.T); % sound velocity
AOA = zeros(N,length(REC.th)); % angle of arrival
AOA_error = AOA; % angle of arrival error
TRAJ.th = zeros(N,length(REC.th)); % trajectory theta/azimuth
% Figures
mainf=figure;
errorf=figure;

for i=1:length(REC.th)
    %% Create shed and room
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
    % adjust axis
    axis equal
    xlim([0 ROOM.xyz(1)])
    ylim([0 ROOM.xyz(2)])

    %% Create line
    TRAJ.y = linspace(-tm,tm,N)' + REC.xyz(2);
    TRAJ.x = ones(N,1)*(REC.xyz(1)+REC.d);
    TRAJ.z = ones(N,1)*REC.xyz(3);
    TRAJ.xyz = [TRAJ.x TRAJ.y TRAJ.z];
    TRAJ.th(:,i) = linspace(45,135,N) + REC.th(i);
    % draw line
    hold on;
    plot(TRAJ.xyz(:,1), TRAJ.xyz(:,2),'b')
    plot(TRAJ.xyz(1,1), TRAJ.xyz(1,2),'b*')
    hold off;

    %% Algorithm
    % read audio
    [AUDIO.y, CH.fs] = audioread(AUDIO.name);

    % correlation range estimation
    MED_N = MED_t * CH.fs;
    CR = round(MED_N);
    % get block sample size
    BLK_N = BLK_t * CH.fs; % get audio samples <=> 100ms

    % == SIM STEREO CONFIG ==
    SS.REC.type = 0; % no surfaces around the microphone
    SS.REC.struct = 'recstruct';
    SS.REC.loc = REC.xyz;
    SS.REC.th = REC.th(i);
    SS.REC.phi = 0;
    SS.REC.mic.dmf = 1;

    SS.ROOM.size = ROOM.xyz;
    SS.ROOM.coeff = [0 0 0];
    SS.ROOM.MR = 0;
    SS.ROOM.R = Room();
    SS.ROOM.R.T = ROOM.T;
    SS.ROOM.R.H = ROOM.H;
    SS.ROOM.R.fs = CH.fs;
    SS.ROOM.R.fl = fl;
    SS.ROOM.R.fh = fh;

    SS.AUDIO.fs = CH.fs;
    SS.AUDIO.fl = fl;
    SS.AUDIO.fh = fh;
    % == EndOfConfig ==

    % predictor loop
    ELAPSED=0;
    for j=1:N
        i1 = 1+(j-1)*BLK_N;
        i2 = j*BLK_N;
        SS.AUDIO.s = AUDIO.y(i1 : i2, 1);
        SPK.loc = TRAJ.xyz(j,:);
        [CH.L, CH.R] = sim_stereo(SS.REC, SS.ROOM, SPK, SS.AUDIO, 1);
        TMR=tic;
        [AOA(j,i),~] = detect_az3(CH, CR, C, REC.DX); %pred v3 (spline interp)
%         [AOA(j,i),~] = detect_az2(CH, CR, C, REC.DX); %pred v2 (3point interp)
    %     AOA(j,i) = detect_az1(CH, CR, C, REC.DX); %pred v1 (non-interp)
        ELAPSED=ELAPSED+toc(TMR);
    end
    ELAPSED = ELAPSED/N;
    fprintf("========== RECEIVER THETA: %dº ==========\n", REC.th(i));
    fprintf("Predictor avg time: %.3f us\n", ELAPSED*1e6);
    fprintf("AOA: min=%.1fº max=%.1fº\n",min(AOA(:,i)),max(AOA(:,i)));
    
    AOA_error(:,i) = abs(AOA(:,i) - TRAJ.th(:,i));
end

%% Adjust labels
figure(mainf)
subplot(2,length(REC.th),1)
xlabel('m','Rotation',0)
ylabel('m','Rotation',0)
subplot(2,length(REC.th),2)
xlabel('m','Rotation',0)
ylabel('m','Rotation',0)

%% Final figures
figure(mainf);
subplot(2,length(REC.th),length(REC.th)+1:2*length(REC.th));
plot(TRAJ.th, AOA,'b')
hold on;
plot(TRAJ.th, TRAJ.th, 'r:')
hold off;
xlabel('AOA Reference (º)')
ylabel('AOA Predictor (º)')

figure(errorf);
subplot(2,1,1)
plot(TRAJ.th, AOA_error, 'b')
subplot(2,1,2)
semilogy(TRAJ.th, AOA_error, 'b')
axh=axes(errorf,'visible','off');
axh.Title.Visible='on';
axh.XLabel.Visible='on';
axh.YLabel.Visible='on';
xlabel(axh,'AOA Reference (º)')
ylabel(axh,'AOA Error (º)')
SIM.TRAJ.th = reshape(TRAJ.th,[],1);
SIM.AOA = reshape(AOA,[],1);
SIM.AOA_error = reshape(AOA_error,[],1);
save('sim_line.mat','-struct','SIM');