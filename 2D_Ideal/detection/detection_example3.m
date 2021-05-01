clear;clc;close all;
addpath ../../simLocUS
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../structures

%% Parameters
% constants
C = 340;
% head porperties
rstruct = 'recstruct.mat';
theta = 0;
phi = 0;
dmf = 1;
diam_head = 0.1*dmf*2;
% room properties
roomc = [5 7 3];
recloc = [1 5 1.7];
spkloc = [2 1 1.7];
% function settings
hidefig = 0;
% interpolation parameters
% range_inter = 2^3;
% inter_points = 2^9;

%% Generate stereo
[y, fs] = audioread('cportugal.wav');
rangemax = round(diam_head + 1e-3*fs);
[yL, yR] = sim_stereo(  rstruct,... 
                        roomc,...
                        recloc,...
                        [theta, phi],...
                        spkloc,...
                        y, fs,...
                        hidefig,...
                        dmf);

%% Draw auxiliary lines
% draw auxiliary horizontal line across receiver
% (to help visualize the azimuth)
if hidefig==0
    point1 = [0 -2 0];
    point2 = [0 2 0];
    point1 = rotate2D(point1, theta);
    point2 = rotate2D(point2, theta);
    lx = [point1(1) point2(1)]+recloc(1);
    ly = [point1(2) point2(2)]+recloc(2);
    lz = [point1(3) point2(3)]+recloc(3);
    hold on;plot3(lx, ly, lz, 'b--');hold off
    % view in 2D
    view(2)
end

%% YAW testing
lr_corr = xcorr(yL,yR, rangemax);
[value, index] = max(abs(lr_corr));
delay_index = (index-1) - rangemax;
delay_t = delay_index*1/fs;
azimuth = acosd(C*delay_t/diam_head)

hold on;
pL = sqrt((recloc(1)-spkloc(1))^2 + (recloc(2)-spkloc(2))^2);
p0 = recloc;
DAz = 90 - azimuth;
aX = p0(1)+pL*(cosd(DAz));
aY = p0(2)-pL*(sind(DAz));
p1 = spkloc;
p2 = [aX aY recloc(3)];

LX1 = [p0(1) p1(1)];
LY1 = [p0(2) p1(2)];
LZ1 = [p0(3) p1(3)];

LX2 = [p0(1) p2(1)];
LY2 = [p0(2) p2(2)];
LZ2 = [p0(3) p2(3)];

plot3(LX1, LY1, LZ1, 'k');
plot3(LX2, LY2, LZ2, 'r');
hold off;