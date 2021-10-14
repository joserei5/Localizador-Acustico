clear force all;
clc;
close force all;
addpath ../functions

%% Audio device info
qc.fs = 48e3;   % sampling frequency
qc.b  = 24;     % audio bits

%% Algorithm Variables
MED_t = 1e-3;                   % maximum estimated delay = 1 ms
MED_N = MED_t * qc.fs;          % maximum estimated delay (in samples)

ROOM.T = 24.3;                  % room temperature (ºC)
C = 20.05*sqrt(273.15 + ROOM.T);% sound velocity (m/s)
REC.d = 29.2e-2;                % receiver distance = 29.5 cm

CH.fs = qc.fs;
CR = round(MED_N);     % correlation range=receiver distance + MED (in samples)

%% Retrieve audio + Prepare blocks
[file,path] = uigetfile('../soundfiles/capture/*.wav');

if isnumeric(file) || isnumeric(path)
    return;
end

snd = join([path,file]);
[y,~] = audioread(snd);

BLK_t = 100e-3;             % Block size = 100 ms
BLK_N = BLK_t * qc.fs;      % 100 blocks of 100ms (10 s of audio)

%% Retrieve experimental results
AOA1 = zeros(100,1);
AOA2 = zeros(100,1);
Delay = zeros(100,1);

for i=1:100
    CH.L = y(1+(i-1)*BLK_N:i*BLK_N,1);
    CH.R = y(1+(i-1)*BLK_N:i*BLK_N,2);
    AOA1(i,1) = detect_az1(CH,CR,C,REC.d);
    [AOA2(i,1), Delay(i,1)]= detect_az2(CH,CR,C,REC.d);
end


%% Non-altered plots
figure;
t=BLK_t:BLK_t:10;
plot(t,AOA1);
hold on;
plot(t,AOA2);
hold off;
xlabel('seconds')
ylabel('AOA (degrees)')
title('AOA(t)')
legend('No Interp.','Interp.')

%% Fix value "spikes"
threshold = 10;
for i=1:100
    if AOA1(i,1) > median(AOA1)+threshold || AOA1(i,1) < median(AOA1)-threshold
        AOA1(i,1) = median(AOA1);
    end
    if AOA2(i,1) > median(AOA2)+threshold || AOA2(i,1) < median(AOA2)-threshold
        AOA2(i,1) = median(AOA2);
    end
end

figure;
t=BLK_t:BLK_t:10;
plot(t,AOA1);
hold on;
plot(t,AOA2);
hold off;
xlabel('seconds')
ylabel('AOA (degrees)')
title('AOA(t)')
legend('No Interp.','Interp.')