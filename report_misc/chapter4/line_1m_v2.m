clear all;clc;close all;
addpath ../../functions
addpath ../../soundfiles/generated
addpath ../../simLocUS

%% Variables
%##Room
ROOM.xyz = [3 3 3]; % room coordinates
ROOM.T = 24; % room temperature
ROOM.H = 75; % room humidity
%##Receiver
REC.xyz = [0.7 ROOM.xyz(2)/2 1.7]; % receiver coordinates
REC.d = 1; % receiver distance in relation to the source
REC.DX = 29.2*1e-2; % receiver: microphone interdistance
REC.th = [-45 45]; % receiver theta
Nexp = length(REC.th);
%##Audio settings
AUDIO.name = {'neg45_line2.wav';'pos45_line2.wav'}; % audio(source) file name
AUDIO.fs = [audioinfo(AUDIO.name{1}).SampleRate;...
            audioinfo(AUDIO.name{2}).SampleRate ];
BLK_t = 100*1e-3; % block time size
BLK_N = BLK_t * AUDIO.fs; % block sample size
N1 = audioinfo(AUDIO.name{1}).TotalSamples / (AUDIO.fs(1)*BLK_t);
N2 = audioinfo(AUDIO.name{2}).TotalSamples / (AUDIO.fs(2)*BLK_t);
N=N1+N2;
Nloop = [N1 N2];
%##Trajectory settings
tm = 1; % trajectory distance from center (90degrees)
%##Algorithm settings
MED_t = 1e-3; % maximum estimated delay
MED_N = MED_t * AUDIO.fs; % correlation range estimation
CR = round(MED_N); % correlation range
C = 20.05*sqrt(273.15+ROOM.T); % sound velocity
AOA = zeros(N,1); % angle of arrival
AOA_error = zeros(N,1); % angle of arrival error
TRAJ.th = zeros(N,1); % trajectory theta/azimuth
%##Figures
mainf=figure;
errorf=figure;

for i=1:Nexp
    % Main loop indexing
    il_1 = N1*(i-1) + 1;
    il_2 = N1 + N2*(i-1);
    il_r = il_1:il_2;
    
    %% Create shed and room
    figure(mainf);
    subplot(2,Nexp,i)
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
    TRAJ.y = linspace(-tm,tm,Nloop(i))' + REC.xyz(2);
    TRAJ.x = ones(Nloop(i),1)*(REC.xyz(1)+REC.d);
    TRAJ.z = ones(Nloop(i),1)*REC.xyz(3);
    TRAJ.xyz = [TRAJ.x TRAJ.y TRAJ.z];
    TRAJ.th(il_r,1) = linspace(45,135,Nloop(i)) + REC.th(i);
    % draw line
    hold on;
    plot(TRAJ.xyz(:,1), TRAJ.xyz(:,2),'b')
    plot(TRAJ.xyz(1,1), TRAJ.xyz(1,2),'b*')
    hold off;

    %% Algorithm
    % read audio
    [AUDIO.y, ~] = audioread(AUDIO.name{i});
    % update algorithm variables
    CH.fs = AUDIO.fs(i);

    % predictor loop
    ELAPSED=0;
    for j=1:Nloop(i)
        % indexing
        i1 = 1+(j-1)*BLK_N(i);
        i2 = j*BLK_N(i);
        iloop = N1*(i-1)+j;
        % get audio channels
        CH.L = AUDIO.y(i1 : i2, 1);
        CH.R = AUDIO.y(i1 : i2, 2);
        
        TMR=tic;
%         [AOA(iloop,1),~] = detect_az3(CH, CR(i), C, REC.DX); %pred v3 (spline interp)
        [AOA(iloop,1),~] = detect_az2(CH, CR(i), C, REC.DX); %pred v2 (3point interp)
%         AOA(iloop,1) = detect_az1(CH, CR(i), C, REC.DX); %pred v1 (non-interp)
        ELAPSED=ELAPSED+toc(TMR);
    end
    ELAPSED = ELAPSED/Nloop(i);
    fprintf("========== RECEIVER THETA: %dº ==========\n", REC.th(i));
    fprintf("Predictor avg time: %.3f us\n", ELAPSED*1e6);
    fprintf("AOA: min=%.1fº max=%.1fº\n",min(AOA(il_r,1)),max(AOA(il_r,1)));
    
    % Calculate AOA errors
    AOA_error(il_r,1) = abs(AOA(il_r,1) - TRAJ.th(il_r,1));
end

%% Final figures
figure(mainf);
subplot(2,Nexp,Nexp+1:2*Nexp);
plot(TRAJ.th, AOA,'b')
hold on;
plot(TRAJ.th, TRAJ.th, 'r:')
hold off;

figure(errorf);
subplot(2,1,1)
plot(TRAJ.th, AOA_error, 'b')
subplot(2,1,2)
semilogy(TRAJ.th, AOA_error, 'b')

SIM.TRAJ.th = TRAJ.th;
SIM.AOA = AOA;
SIM.AOA_error = AOA_error;
save('sim_line_v2.mat','-struct','SIM');