%% Generate stereo
clear;clc;close all;
[y, fs] = audioread('cportugal.wav');

theta = -0;
phi = 0;
recloc = [2 5 1.7];
spkloc = [4 5 1.7];
hidefig = 0;
[yL, yR] = sim_stereo([5 7 3], recloc, [theta, phi], spkloc, y, fs, hidefig);

diam_head = 0.5;

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
    hold on;plot3(lx, ly, lz, 'b', 'LineWidth', 1.5);hold off
    % view in 2D
    view(2)
end

%% Time unit test simulation
%{
    1000 runs; each run = mean of 1000 experiences
    delayseq(yL,40); = 0.002133609477200s
    dot(yL,yR) = 1.873728969999999e-05s
    [zeros(40,1); yL(1+40:end)] ~ 0.0014s
    ressss = circshift(yL,40); \ 
    ressss(1:40) = 0;          / ~0.0007s
%}

%% Plot channels
t = (0:length(y)-1)*1/fs;
rangemax = round(1e-3*fs);

lr_corr = xcorr(yL,yR, rangemax);

figure;
plot(-rangemax:rangemax, abs(lr_corr))

[value, index] = max(abs(lr_corr));
delay_index = (index-1) - rangemax;
delay_t = delay_index*1/fs;
azimuth = acosd(340*delay_t/diam_head);

fprintf("Delay of %.1f us | ", delay_t*1e6)
fprintf("Azimuth = %.1fº\n", 180-azimuth)


%% Lagrange interpolation
range_inter = 2^2;
inter_points = 2^9;

x_lg = index-range_inter:index+range_inter;
y_lg = lr_corr(x_lg);

x_interp = linspace(min(x_lg), max(x_lg), inter_points);
y_interp = interp1(x_lg, y_lg, x_interp, 'spline');

figure;
subplot(211)
plot(x_lg, y_lg);
subplot(212)
plot(x_interp, y_interp);

[value, index] = max(y_interp);
delay_index = x_interp(index-1) - rangemax;
delay_t = delay_index*1/fs;
azimuth = acosd(340*delay_t/diam_head);

fprintf("Delay of %.1f us | ", delay_t*1e6)
fprintf("Azimuth = %.1fº\n", 180-azimuth)


%% Simple interpolation
% t = (0:length(y)-1)*1/fs;
% rangemax = round(1e-3*fs);
% 
% lr_corr = xcorr(yL,yR, rangemax);
% [value, index] = max(abs(lr_corr));
% p = index-1:index+1;
% 
% c = lr_corr(p(2));
% a = lr_corr(p(1)) + lr_corr(p(3)) - 2*c;
% b = a + c - lr_corr(p(1));
% 
% x = -b/(2*a);
% disp(value)
% disp(value-x)
% 
% r = roots([a b c-x])

%% Stereo sound
% y = [yL'; yR'];
% soundsc(y,fs)

%% Compare interpolation results
theta_range = -90:0.5:90;
N = length(theta_range);
azimuth_l = zeros(N, 4);

rangemax = round(1e-3*fs);

range_inter = 2^3;
inter_points = 2^9;

for i=1:N
    [yL, yR] = sim_stereo([5 7 3], recloc, [theta_range(i), 0], spkloc, y, fs, 1);
    
    lr_corr = xcorr(yL,yR, rangemax);
    [value, index] = max(abs(lr_corr));
    delay_index = index - rangemax;
    delay_t = delay_index*1/fs;
    azimuth_l(i,1) = acosd(340*delay_t/diam_head);
    
%     lim_inf = index - range_inter;
%     lim_sup = index + range_inter;
%     
%     if lim_inf <= 0
%         lim_inf = 1;
%     end
%     
%     if lim_sup > length(lr_corr)
%         lim_sup = length(lr_corr);
%     end
%     
%     x_lg = lim_inf:lim_sup;
%     y_lg = lr_corr(x_lg);
%     x_interp = linspace(min(x_lg), max(x_lg), inter_points);
%     
%     y_interp = interp1(x_lg, y_lg, x_interp, 'spline');
%     [value, index] = max(y_interp);
%     delay_index = x_interp(index) - rangemax;
%     delay_t = delay_index*1/fs;
%     azimuth_l(i,2) = acosd(340*delay_t/diam_head);
%     
%     y_interp = interp1(x_lg, y_lg, x_interp, 'pchip');
%     [value, index] = max(y_interp);
%     delay_index = x_interp(index) - rangemax;
%     delay_t = delay_index*1/fs;
%     azimuth_l(i,3) = acosd(340*delay_t/diam_head);
%     
%     y_interp = interp1(x_lg, y_lg, x_interp, 'makima');
%     [value, index] = max(y_interp);
%     delay_index = x_interp(index) - rangemax;
%     delay_t = delay_index*1/fs;
%     azimuth_l(i,4) = acosd(340*delay_t/diam_head);
%     
end

figure;
subplot(221)
    plot(theta_range, theta_range+90, 'k')
    hold on
    plot(theta_range, azimuth_l(:,1), 'b')
    hold off
    title('No interpolation')
    xlabel('theta (degrees)')
    ylabel('predicted azimuth (degrees)')
% subplot(222)
%     plot(theta_range, theta_range+90, 'k')
%     hold on
%     plot(theta_range, azimuth_l(:,2), 'r')
%     hold off
%     title('''spline'' interpolation')
%     xlabel('theta (degrees)')
%     ylabel('predicted azimuth (degrees)')
% subplot(223)
%     plot(theta_range, theta_range+90, 'k')
%     hold on
%     plot(theta_range, azimuth_l(:,3), 'r')
%     hold off
%     title('''pchip'' interpolation')
%     xlabel('theta (degrees)')
%     ylabel('predicted azimuth (degrees)')
% subplot(224)
%     plot(theta_range, theta_range+90, 'k')
%     hold on
%     plot(theta_range, azimuth_l(:,4), 'r')
%     hold off
%     title('''makima'' interpolation')
%     xlabel('theta (degrees)')
%     ylabel('predicted azimuth (degrees)')