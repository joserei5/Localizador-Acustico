%% Initial cmds
clc
clear
close all

%% Parameters
Nx = 100;              % no. points for Microphone Spacing (d_x)
Nd =  80;              % no. points for distance (d)
Ne = 2000;             % no. points for Epsilon range (\epsilon)

epsilon.min = -500e-6;
epsilon.max = +500e-6;

room.mic.dX.min = 0;
room.mic.dX.max = 2;
room.temp = 20;

src.d.min = 0;
src.d.max = 7;
% src.theta = 0;

%% Execute
iterParam1(room, src, epsilon, Nx, Nd, Ne)
