clear; clc; close all;

%% Init variables
MR = 2; % max reflections (order)
fs = 192e3; % sampling frequency (>= 4xfh)
fl = 100;  % lower frequency bound
fh = 60e3; % upper frequency bound (<=80 kHz due to @KemoL10_TF)

%% II) Wall p1
AW = addWall([],[0 0 0], [7 0 0], [7 0 7], [0 0 7], 0.7);

%% II) Acoustic sensors
close all;
AM = addMic([], [9 2.5 3.5], @omni, [1 0 0]);
AM = addMic(AM, [9 2.5 3.5], @omni, [1 0 0]);

makeFile('exampleobj', AW, AM, MR, fs, fl, fh); 
displayRoom('exampleobj','HideVS');

%% II) Complete room
close all;
AW = addDivision([7 7 7], [0 0 0.5]);
makeFile('exampleobj', AW, AM, MR, fs, fl, fh); 
displayRoom('exampleobj','HideVS');

%% II) Microphone p1
close all;
AM = addMic([], [0 0 0], @omni, [1 0 0]);
AM = addMic(AM, [0 1 0], @omni, [1 0 0]);
makeFile('exampleobj', [], AM, MR, fs, fl, fh); 
displayRoom('exampleobj','HideVS');

%% II) Virtual Sources
close all;
AW = addDivision([7 7 7], [0.7 0.8 0.5]);

AM = addMic([], [2 3 3]);
AM = addMic(AM, [2 4 3]);
makeFile('exampleobj', AW, AM, MR, fs, fl, fh); 
displayRoom('exampleobj');
displayRoom('exampleobj','HideVS');

%% III) Wall in the middle
close all;
AW = addDivision([7 7 7], [0 0 0.5]);

AM = addMic([], [2 3.4 3]);
AM = addMic(AM, [2 3.6 3]);

% wall middle
AW = addWall(AW,[1.9 3.5 2.9], [2.1 3.5 2.9], [2.1 3.5 3.1], [1.9 3.5 3.1], 0);

makeFile('exampleobj', AW, AM, MR, fs, fl, fh); 
displayRoom('exampleobj','HideVS');

S = addSpk([6 6 3]);

R = Room(); % Create an Room with typical characteristics
R.T = 25;  % temperatura ºC
R.H = 30;  % humidade %
R.P = 1.01;% pressure atm

I = impR('exampleobj', S, R);


% for convenience lets consider the first mic the left channel and the
% second mic. the right channel 
t = (0:length(I)-1)/fs;
figure
plot(1000*t,abs(I(1,:)),1000*t,abs(I(2,:)))
title('Impulse Response')
ylabel('Amplitude')
xlabel('Time (ms)')
legend('Left Channel','Right Channel')

%% III) Prism-like object - fixed positions
close all;
% head
AW = [];
[AM, AW] = addHumanHead(AW, 2, 3.5, 3, '+xx', 0);
makeFile('exampleobj', AW, AM, MR, fs, fl, fh); 
displayRoom('exampleobj','HideVS');

AW = addDivision([7 7 7], [0 0 0.5]);
[AM, AW] = addHumanHead(AW, 2, 3.5, 3, '+xx', 0);
makeFile('exampleobj', AW, AM, MR, fs, fl, fh); 
displayRoom('exampleobj','HideVS');

S = addSpk([6 6 3]);

R = Room(); % Create an Room with typical characteristics
R.T = 25;  % temperatura ºC
R.H = 30;  % humidade %
R.P = 1.01;% pressure atm

I = impR('exampleobj', S, R);


% for convenience lets consider the first mic the left channel and the
% second mic. the right channel 
t = (0:length(I)-1)/fs;
figure
plot(1000*t,abs(I(1,:)),'LineWidth', 1.5)
hold on
plot(1000*t,abs(I(2,:)))
hold off
title('Impulse Response')
ylabel('Amplitude')
xlabel('Time (ms)')
legend('Left Channel','Right Channel')