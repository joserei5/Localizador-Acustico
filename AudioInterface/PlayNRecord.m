clear force all;
clc;
close force all;
addpath ../soundfiles/generic/

%% Get drivers
info = audiodevinfo;
driversIN = cell(length(info.input),1);
for i=1:length(info.input)
    driversIN{i,1} = join([string(info.input(i).ID), ' - ', info.input(i).Name]);
end
driversIN

driversOUT = cell(length(info.output),1);
for i=1:length(info.output)
    driversOUT{i,1} = join([string(info.output(i).ID), ' - ', info.output(i).Name]);
end
driversOUT

%% Prepare audio devices
qc.fs = 48e3;   % sampling frequency
qc.b = 24;      % audio bits
qc.dev.r = 1;   % 1-2 (QUAD-CAPTURE) (Windows DirectSound) -- input
qc.dev.p = 11;   % 1-2 (QUAD-CAPTURE) (Windows DirectSound) -- output
% recorder
recObj = audiorecorder(qc.fs, qc.b, 2, qc.dev.r)
% player -- frequency beep
t = 0:1/qc.fs:0.25;
fbeep = 0.1*sin(2*pi*250*t);
beepObj = audioplayer(fbeep, qc.fs, qc.b, qc.dev.p)
% player -- drone noise
% [drone, dfs] = audioread('drone1.mp3');
% droneObj = audioplayer(drone, dfs, qc.b, qc.dev.p)
% player -- mosquito noise
[mosquito, dfs] = audioread('mosquito2.wav');
mosquitoObj = audioplayer(mosquito, dfs, qc.b, qc.dev.p)


%% Record Audio

% WARNING LOOP
LAPS = 0;
BEEPS = 10;
while LAPS < BEEPS
   % reset timer
   TMR = tic;
    
   % play beep player
   play(beepObj);
   
   % wait approx. 1 second
   while toc(TMR) < 0.95
   end
   
   % force stop beep player
   stop(beepObj);
  
   % LAP
   LAPS = LAPS + 1;
end
clear TMR;

% SEND DESIRED SAMPLE TO OUTPUT
% play(droneObj);
play(mosquitoObj);

% START RECORDING
recordblocking(recObj, 10);

% FORCE STOP OUTPUT
% stop(droneObj);
stop(mosquitoObj);

% FINAL BEEP
play(beepObj);

%% Graphs
y = getaudiodata(recObj);
t=0:1/qc.fs:(length(y)-1)/qc.fs;
t=t';
figure;plot(t,y)

%% Listen to audio
% sound(10*y,qc.fs)

%% Save audio
fname = join(['../soundfiles/capture/','recording',datestr(now,'ddmmyy_HHMMSS'),'.wav']);
audiowrite(fname,y,qc.fs);