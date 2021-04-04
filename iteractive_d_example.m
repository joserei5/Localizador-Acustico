clear;clc;close all;
addpath('wav')

room.size = [5 7 3];
room.rec.loc = [2 5 1.7];
room.spk.loc = [4 4 1.7];
room.temp = 20;

rec.struct = 'recstruct.mat';
rec.mic.dmf = 1;
rec.mic.dist = 0.1*rec.mic.dmf*2;
rec.yaw = 0;
rec.pitch = 0;

distance = 1.5;
points = 1000;

iterCDetection(room, rec, distance, points)