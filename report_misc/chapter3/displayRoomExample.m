clear all; close all; clc;
addpath ../../simLocUS/
addpath ../../functions/

%% Walls Array
cmatrix = [5 5 3]; % rectangular shape room -> size coordinates matrix: [x y z]
rcmatrix = [0.7 0.8 0.5]; % reflection coefficients matrix: [wall ceilling floor]
AW = addDivision(cmatrix, rcmatrix);


%% Receiver Array
% (only using a single microhpone for now)
AM = addMic([], [2.5 0.5 1.0], @omni, [0 1 0]);

%% Room Representation File: No reflections
MR = 0; % 
fs = 4*48e3; % sampling frequency (>= 4xfh)
fl = 100;  % lower frequency bound
fh = 20e3; % upper frequency bound (<=80 kHz due to @KemoL10_TF)

makeFile('testfile', AW, AM, MR, fs, fl, fh);

%% Display the Room: No reflections
displayRoom('testfile');
set(gcf,'Renderer','painters');
% displayRoom('testfile','HideVS');

%% Room Representation File: With MR=1
MR = 1; %
makeFile('testfile', AW, AM, MR, fs, fl, fh);

%% Display the Room: With MR=1
displayRoom('testfile');
set(gca,'Color',[0.8 0.8 0.8])
set(gcf,'Renderer','painters');
% displayRoom('testfile','HideVS');

%% Room Representation File: With MR=2
MR = 2; %
makeFile('testfile', AW, AM, MR, fs, fl, fh);

%% Display the Room: With MR=2
displayRoom('testfile');
set(gca,'Color',[0.8 0.8 0.8])
set(gcf,'Renderer','painters');
% displayRoom('testfile','HideVS');

%% Room Representation File: With MR=4
MR = 4; %
makeFile('testfile', AW, AM, MR, fs, fl, fh);

%% Display the Room: With MR=4
displayRoom('testfile');
set(gca,'Color',[0.8 0.8 0.8])
set(gcf,'Renderer','painters');
% displayRoom('testfile','HideVS');