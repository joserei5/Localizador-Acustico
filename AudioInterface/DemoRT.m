fclose('all');
clear force all;
clc;

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
qc.dev.p = 7;   % 1-2 (QUAD-CAPTURE) (Windows DirectSound) -- output

% recorder
recObj = audiorecorder(qc.fs, qc.b, 2, qc.dev.r);
set(recObj,'TimerPeriod',100e-3,'TimerFcn',{@fff,recObj});

stop(recObj)
record(recObj)

%% Functions
function fff(~,~,recObj)
    toc
    timefff;
%     stop(recObj);
%     Y = getaudiodata(recObj);
%     record(recObj);
end

function timefff(~,~)
    tic
end