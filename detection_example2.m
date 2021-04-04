clear;clc;close all;
addpath('wav')

%% Parameters
% constants
C = 340;
% head porperties
theta = 0;
phi = 0;
dmf = 1;
diam_head = 0.1*dmf*2;
% room properties
recloc = [2 5 1.7];
spkloc = [4 5 1.7];
% function settings
hidefig = 0;
% interpolation parameters
% range_inter = 2^3;
% inter_points = 2^9;

%% Generate stereo
[y, fs] = audioread('cportugal.wav');
rangemax = round(1e-3*fs);
[yL, yR] = sim_stereo(  'recstruct.mat',...
                        [5 7 3],...
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
% theta range array
theta_range = -90:0.5:90;
N = length(theta_range);
% predicted theta array
azimuth_l = zeros(N, 4);
% distance range array
d_range = 1:10;
K = length(d_range);

for k=1:K
    diam_head = 0.1*d_range(k)*2;
    for i=1:N
        [yL, yR] = sim_stereo('recstruct.mat',...
                            [5 7 3],...
                            recloc,...
                            [theta_range(i), phi],...
                            spkloc,...
                            y, fs,...
                            1,...
                            d_range(k));

        lr_corr = xcorr(yL,yR, rangemax);
        [value, index] = max(abs(lr_corr));
        delay_index = (index-1) - rangemax;
        delay_t = delay_index*1/fs;
        
        if abs(C*delay_t/diam_head) > 1
            disp(C*delay_t/diam_head);
        end
        azimuth_l(i,1) = acosd(C*delay_t/diam_head);
    end

    fprintf('%.0d ', k);
    figure;
        plot(theta_range, theta_range+90, 'k')
        hold on
        plot(theta_range, real(azimuth_l(:,1)), 'b')
        hold off
        addinfo = (['microphone distance = ', num2str(diam_head), 'm']);
        title('No interpolation', addinfo)
        xlabel('theta (degrees)')
        ylabel('predicted azimuth (degrees)')
end



% ttt = -90:15:90;
% for i=1:length(ttt)
%    [yL, yR] = sim_stereo('recstruct.mat',...
%                             [5 7 3],...
%                             recloc,...
%                             [ttt(i), phi],...
%                             spkloc,...
%                             y, fs,...
%                             0,...
%                             10);
%     point1 = [0 -2 0];
%     point2 = [0 2 0];
%     point1 = rotate2D(point1, theta);
%     point2 = rotate2D(point2, theta);
%     lx = [point1(1) point2(1)]+recloc(1);
%     ly = [point1(2) point2(2)]+recloc(2);
%     lz = [point1(3) point2(3)]+recloc(3);
%     hold on;plot3(lx, ly, lz, 'b', 'LineWidth', 1.5);hold off
%     % view in 2D
%     view(2)
% end