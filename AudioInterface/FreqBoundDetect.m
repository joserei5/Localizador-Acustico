clear; clc; close all;
addpath ../functions/

%% Variables
% show sound wave amplitude (time)
F_SWAVEAMP = 1;
% show frequency plot
F_FREQPLOT = 1;
% show spectrograms
F_SPECTROGRAM = 1;
% show filter visualisation tool
F_FVTOOL = 1;

% real-time OR local sound file
% scenario (flag)
% 0 = offline; 1 = real-time;
F_offlineRT = 0;

% music file dir
% M.dir = '../soundfiles/generic/mosquito.wav';
M.dir = '../soundfiles/capture/extensive/1m/deg90.wav';

% block size
% 100 ms
B.t = 100*1e-3;

% type of filter (flag)
% 0 = lowpass;
% 1 = highpass;
% 2 = bandpass;
F_TYPE = 2;
% filter order
F.order = 200;
% filter frequency bandpass cut-off values
F.bp.co1 = 190;
F.bp.co2 = 1e3;
% filter frequency highpass cut-off value
F.hp.co = 200;
% highpass filter stopband attenuation
F.hp.spAtt = 80;
% filter frequency lowpass cut-off value
F.lp.co = 10e3;
% filter sample rate
F.fs = 48e3;

% hamming window size
Nspec = 256;

% algorithm Variables:
% maximum estimated delay = 1 ms
MED_t = 1e-3;
% room temperature (ºC)
ROOM.T = 24.5;
% sound velocity (m/s)
C = 20.05*sqrt(273.15 + ROOM.T);
% receiver distance = 29.5 cm
REC.d = 29.2e-2;

%% Process Variables
% read sound file
[M.y, M.fs] = audioread(M.dir);

% get respective reference
% AOA value printed on the file name
M.ref=split(M.dir,'/');
M.ref=split(M.ref(end),'.');
M.ref=regexp(M.ref,'\d*','Match');
M.ref=str2num(cell2mat(M.ref{1}));

% algorithm information:
% channel sampling frequency
CH.fs = M.fs;
% maximum estimated delay (samples)
MED_N = MED_t * CH.fs;
% correlation range = receiver distance + MED (samples)
CR = round(MED_N);

% absolute value of fft of sound file values
Ymono = M.y(:,1);
YMONO = abs(fft(Ymono));

% filter type and object
switch F_TYPE
    case 0
        F.type = 'lowpassfir';
        bpFilt = designfilt(    F.type,'FilterOrder',F.order, ...
                                'CutoffFrequency',F.lp.co,...
                                'SampleRate',F.fs);
    case 1
        F.type = 'highpassfir';
        bpFilt = designfilt(    F.type,'FilterOrder',F.order, ...
                                'CutoffFrequency',F.hp.co,...
                                'StopbandAttenuation', F.hp.spAtt,...
                                'SampleRate',F.fs);
    case 2
        F.type = 'bandpassfir';
        bpFilt = designfilt(    F.type,'FilterOrder',F.order, ...
                                'CutoffFrequency1',F.bp.co1,...
                                'CutoffFrequency2',F.bp.co2, ...
                                'SampleRate',F.fs);
end

% filtered sound file values                
Ymono_filt = fftfilt(bpFilt, Ymono);
% absolute value of fft of filtered sound file values
YMONO_FILT = abs(fft(Ymono_filt));

% block sample size
B.N = B.t * M.fs;
% total blocks
B.no = length(Ymono)/B.N;

% AOA list
% experimental
AOA.exp.f2 = zeros(B.no, 1);
AOA.exp.f2f = zeros(B.no, 1);
AOA.exp.f3 = zeros(B.no, 1);
AOA.exp.f3f = zeros(B.no, 1);

% hamming window
wspec = hamming(Nspec);
% Number of overlap samples
Noverlap = Nspec/2;

%% Example Figures
if F_FVTOOL
    fvtool(bpFilt)
end

if F_SWAVEAMP
    Figs.h1 = figure();
    plot(Ymono)
    hold on;
    plot(Ymono_filt)
    hold off;
    xlabel('t (s)');xlabel('A');
end

if F_FREQPLOT
    Figs.h2 = figure();
    subplot(2,1,1)
    stem(YMONO)
    subplot(2,1,2)
    stem(YMONO_FILT)
    xlabel('f (Hz)');ylabel('N');
end

if F_SPECTROGRAM
    Figs.h3 = figure();
%     [Ssteady, fspec, tspec] = spectrogram(M.y, wspec, Noverlap, Nspec, M.fs);
    spectrogram(Ymono, wspec, Noverlap, Nspec, M.fs, 'yaxis');
    Figs.h4 = figure();
    spectrogram(Ymono_filt, wspec, Noverlap, Nspec, M.fs, 'yaxis');
end
drawnow

%% Batch processing of saved soundfile
if ~F_offlineRT
    for cnt = 1:B.no
        idx1 = (B.N-1)*cnt + 1;
        idx2 = B.N*cnt;
        CH.L = M.y(idx1:idx2,1);
        CH.R = M.y(idx1:idx2,2);
        
        [AOA.exp.f2(cnt,1), ~] = detect_az2(CH, CR, C, REC.d);
        [AOA.exp.f3(cnt,1), ~] = detect_az3(CH, CR, C, REC.d);
        [AOA.exp.f2f(cnt,1), ~] = detect_az2_filtered(CH, CR, C, REC.d, bpFilt);
        [AOA.exp.f3f(cnt,1), ~] = detect_az3_filtered(CH, CR, C, REC.d, bpFilt);
    end
end


%% Final Figures

Figs.h5 = figure();
plot(AOA.exp.f2,'b')
hold on;
plot(AOA.exp.f3,'b*')
plot(AOA.exp.f2f,'r')
plot(AOA.exp.f3f,'r*')
yline(M.ref,'k--')
hold off;
ylabel('Degrees')
tstring = sprintf('Experimental results for AOA=%.0d degrees', M.ref);
title(tstring)
legend( 'detect\_az2','detect\_az3',...
        'detect\_az2\_filtered','detect\_az3\_filtered',...
        'AOA reference')