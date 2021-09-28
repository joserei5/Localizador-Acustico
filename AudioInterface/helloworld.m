clear force all;
clc;
close force all;

%% Get drivers (INPUT)
info = audiodevinfo;
driversIN = cell(length(info.input),1);
for i=1:length(info.input)
    driversIN{i,1} = join([string(info.input(i).ID), ' - ', info.input(i).Name]);
end
driversIN

%% Prepare recorder
qc.f = 48e3;
qc.b = 24;
qc.dev = 1; % 1-2 (QUAD-CAPTURE) (Windows DirectSound)
recObj = audiorecorder(qc.f, qc.b, 2, qc.dev)

timer=tic;
regtimer=0.95;
i=1;
while toc(timer)<3.0
   fprintf("%d ",i);
   i=i+1;
   while toc(timer) < regtimer
   end
   regtimer = regtimer + toc(timer);
end
clear timer;

disp('Start speaking.')
recordblocking(recObj, 10);
disp('End of Recording.');

%% Graphs
y = getaudiodata(recObj);
t=0:1/qc.f:(length(y)-1)/qc.f;
t=t';
figure;plot(t,y)

%% Save data
REC.obj = recObj;
REC.t=t;
REC.y=y;
REC.qc = qc;
save('recdata','-struct','REC');
audiowrite('recording.wav', REC.y, REC.qc.f, 'BitsPerSample', REC.qc.b);

%% Get drivers (INPUT)
info = audiodevinfo;
driversOUT = cell(length(info.output),1);
for i=1:length(info.output)
    driversOUT{i,1} = join([string(info.output(i).ID), ' - ', info.output(i).Name]);
end
driversOUT

%% Load data
clear force all;
clc;
close force all;
REC = load('recdata.mat');
figure; plot(REC.t, REC.y)


player = audioplayer(REC.y,REC.qc.f,REC.qc.b,7);
play(player)

%% Algorithm Variables
MED_t = 1e-3;                   % maximum estimated delay = 1 ms
MED_N = MED_t * REC.qc.f;       % maximum estimated delay (in samples)

ROOM.T = 27;                    % room temperature (ºC)
C = 20.05*sqrt(273.15 + ROOM.T);% sound velocity (m/s)
REC.d = 20e-3;                  % receiver distance = 20 cm
REC.dN = (REC.d/C)*REC.qc.f;    % receiver distance (in samples)

CH.fs = REC.qc.f;

CR = round(REC.dN + MED_N);     % correlation range=receiver distance + MED (in samples)

%% AudioWriter
snd = 'recording.wav';      % recording directory
BLK_t = 100e-3;             % Blok size = 100 ms
BLK_N = BLK_t * REC.qc.f;   % Block size (in samples)

fileReader = dsp.AudioFileReader(snd);
fileReader.SamplesPerFrame = BLK_N;

fileInfo = audioinfo(snd)
deviceWriter = audioDeviceWriter('SampleRate',fileInfo.SampleRate);
bufferLatency = fileReader.SamplesPerFrame/deviceWriter.SampleRate

setup(deviceWriter,zeros(fileReader.SamplesPerFrame,fileInfo.NumChannels))
info(deviceWriter)

%% Loop FFT preparation
figure;
d.f=[1];
d.fft.L=[1];
d.fft.R=[1];
subplot(211);SP1=stem(d.f, d.fft.L);
subplot(212);SP2=stem(d.f, d.fft.R);
SP1.XDataSource = 'd.f';
SP1.YDataSource = 'd.fft.L';
SP2.XDataSource = 'd.f';
SP2.YDataSource = 'd.fft.R';


%% Detection Loop
addpath ../functions
AOA = zeros(fileInfo.TotalSamples/BLK_N,1);

i=1;
while ~isDone(fileReader)
    audioData = fileReader();
    
    CH.L = audioData(:,1);
    CH.R = audioData(:,2);
    
    AOA(i,1) = detect_az1(CH,CR,C,REC.d);
    i=i+1;
    
    d.f = (1:length(CH.L))*22e3/length(CH.L);
    d.fft.L = abs(fft(CH.L));
    d.fft.R = abs(fft(CH.R));
    
    refreshdata;
    drawnow;
    
    deviceWriter(audioData);
end

release(fileReader)
release(deviceWriter)