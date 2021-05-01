close all;
clear all;
clc;
% Copyright 2020. 
% Created by Daniel Filipe Albuquerque.
% Minor additions by José Reis

%% LOADING/DEFINING RECEIVER STRUCTURE
load recstruct.mat
% rec.s.pos = wall coordinates of receiver
% rec.s.rc = wall reflection coefficients
% rec.mic.pos = microphone position coordinates
% rec.mic.dir = microphone direction coordinates
% rec.mic.tf = microphone transfer function (only 1x1 dimension)

% trying to implement circulr piston TF
rec.mic.tf = @omni;

%% STAGE ONE
AW = addDivision([5 7 3],[0.7 0.8 0.5]);
[AM, AW] = addReceiver(AW, rec, [1.0 5.5 1.7], [0 0]);

MR = 0; % max reflections (order)
fs = 48e3; % sampling frequency (>= 4xfh)
fl = 100;  % lower frequency bound
fh = 20e3; % upper frequency bound (<=80 kHz due to @KemoL10_TF)

makeFile('headwall', AW, AM, MR, fs, fl, fh);

% AM = [];
% addMic(AM, );
% addMic(AM, );
%% DISPLAY THE ROOM
%displayRoom('headwall');    
displayRoom('headwall','HideVS');
S = addSpk( [4.0 1.0 1.7]);
hold on; plot3(4.0, 1.0, 1.7, 'k*', 'LineWidth', 3); hold off;

%% STAGE 2


R = Room();

R.T = 25;  % temperatura ºC
R.H = 30;  % humidade %
R.P = 1.01;% pressure atm

I = impR('headwall', S, R); % Compute the impulse response for each microphone
t = (0:length(I)-1)/fs;

figure
plot(1000*t,abs(I(1,:)),1000*t,abs(I(2,:)))
title('Impulse Response')
ylabel('Amplitude')
xlabel('Time (ms)')
legend('Left Channel','Right Channel')

%% DISPLAY AUDIO CHANNELS + AUDIO OUTPUT
% Input signal:
%[x, f_x] = audioread('robotvoice.mp3');
[x, f_x] = audioread('cportugal.wav');
%[x, f_x] = audioread('mosquito.wav');
x=x/100000;
% Output signal left channel:
yL = fftfilt(I(1,:),x);
% Output signal right channel:
yR = fftfilt(I(2,:),x);
t = (0:length(x)-1)*1/f_x;

figure
plot(1000*t,yL,1000*t,yR)
title('Signals')
xlim([1000 1010])
ylabel('Amplitude')
xlabel('Time (ms)')
legend('Left Channel','Right Channel')

y = [yL'; yR'];
soundsc(y,f_x)
