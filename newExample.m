close all;
clear all;
clc;
%
%  Copyright 2020.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 1.1$    $Date: 2020/10/30$
%
%% STAGE ONE
%
% In this stage we only define the room to simulate and place
% the microphones inside it. At the end of this stage we save all the
% data for next stage. We only need to do this once.
%
% This example considers a room with four walls, one floor and one ceil.
% The reflection coefficient for the walls is 0.7, for the ceil is 0.8 and
% for the floor is 0.5. 
%
% To add a wall to the simulator we use the function addWall. In the first 
% call we use this function like that:
AW = addWall(   [], [0 0 0], [0 0 3], [5 0 3], [5 0 0], 0.7);
%|              |      |________|________|________|       |
% This is the   Use this     |                            This is the 
% array of      only on      This is the coordinates      reflection
% surfaces      first call   of the vertices of the wall  coefficient
%
% We add the other tree walls:
AW = addWall(AW, [0 0 0], [0 0 3], [0 5 3], [0 5 0], 0.7);
%                 |
%                 Note: if is not the first call you must give to the
%                 function the array that you want to add the surface (AW
%                 in this case)
AW = addWall(AW, [5 0 0], [5 0 3], [5 5 3], [5 5 0], 0.7);
AW = addWall(AW, [0 5 0], [0 5 3], [5 5 3], [5 5 0], 0.7);

% We add the floor and ceil with the respective reflection coefficients
AW = addWall(AW, [0 0 3], [5 0 3], [5 5 3], [0 5 3], 0.8);
AW = addWall(AW, [0 0 0], [5 0 0], [5 5 0], [0 5 0], 0.5);


% To add a microphone to the simulator we use the function addMic. In the
% first call we use this function like that if you have beam function 
% (please see the Speaker part):
% AM = addMic(     [], [2.5 1.0 1.0]          , hb, [1 0 0]);
% |                 |           |                |    |_|_|  
% This is the       Use this    The coordinates  |        |        
% array of          only on     of the source    |        propagation direction
% sources           irst call                    | 
%                                                transfer function
%
% for an omnidirectional microphone you can use:
AM = addMic(     [], [1.0 4.6 1.0]);
% |               |           |                 
% This is the     Use this    The coordinates         
% array of        only on     of the source   
% mics            first call                 
%  

AM = addMic(AM, [1.0 4.8 1.0]);
%            |
%            Note, if is not the first call you must give to the
%            function the array what you want to add the microphone.

% we define that the maximum reflection coefficient that the simulator
% takes into account - MR.
MR = 2; % 
fs = 192e3; % sampling frequency (>= 4xfh)
fl = 100;  % lower frequency bound
fh = 60e3; % upper frequency bound (<=80 kHz due to @KemoL10_TF)

% Now, for finish the stage one, it is necessary to create the file. In
% this example we call the file "testfile". And for this we use the
% function:
makeFile('testfile', AW, AM, MR, fs, fl, fh);
%


%% DISPLAY THE ROOM
%
% You can view in 3D the room, the real sources, and all the virtual
% sources with reflection coefficient up to maxR, for that you can use the
% function displayRoom:
displayRoom('testfile');
displayRoom('testfile','HideVS');
%                        |
%                        This option hide the virtual sources


%% STAGE 2

% To add a Speaker to the simulator we use the function addSpk. 
% It is only possible to add just one Speaker.

% Hear we define the beam function of the Speaker.
%
% It must be a function of type [h,p] = func(theta,phi,fs,fL,fH); where h is
% the impulse response of the wall and p the position of the first impulse
% response sample; theta is the azimuth; phi is the elevation, fs  is the
% sampling frequency and fL and fH are the lower and upper frequency bounds
% for band optimization
%
% hb = @KemoL10_TF; % This is the beam function of a speaker that we use in
% our LAB 

% To add a Speaker to the simulator we use the function addSpk like that:
% S = addSpk( [4.0 2.5 1.0]          , hb, [-1 0 0]);
%                     |                |     |_|_|  
%                     The coordinates  |         |        
%                     of the source    |         propagation direction
%                                      | 
%                                      transfer function

% for an omnidirectional Speaker you can use:
S = addSpk( [1.0 1.0 0.5]);

R = Room(); % Create an Room with typical characteristics
% Default:
% T = 22ºC; P = 1 atm; H = 50%; 

R.T = 25;  % temperatura ºC
R.H = 30;  % humidade %
R.P = 1.01;% pressure atm

I = impR('testfile', S, R); % Compute the impulse response for each microphone
% The impulse reponses are obtained in the order that the microphones were
% added

% for convenience lets consider the first mic the left channel and the
% second mic. the right channel 

t = (0:length(I)-1)/fs;

figure
plot(1000*t,abs(I(1,:)),1000*t,abs(I(2,:)))
title('Impulse Response')
ylabel('Amplitude')
xlabel('Time (ms)')
legend('Left Channel','Right Channel')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plot(1000*t,yR,1000*t,yL)
title('Signals')
xlim([1000 1010])
ylabel('Amplitude')
xlabel('Time (ms)')
legend('Right Channel','Left Channel')

y = [yL'; yR'];
soundsc(y,f_x)
