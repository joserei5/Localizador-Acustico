clear;close all; clc;

%% INITIALIZATION OF VALUES
MR = 2; % max reflections (order)
fs = 192e3; % sampling frequency (>= 4xfh)
fl = 100;  % lower frequency bound
fh = 60e3; % upper frequency bound (<=80 kHz due to @KemoL10_TF)

%% CIRCULAR PISTON
% theta array size
N=2000;
% a-prior known size of frequency range
H_size = 1000;
% frequency range, generated again (1xH_size)
freq_sweep = linspace(0,fs/2,H_size);
% theta range (1xN)
theta_sweep = linspace(-pi, pi, N);
% frequency response array initialization (NxH_size)
theta_h = zeros(N,H_size);

% retrieve values for the circular piston
for i=1:N
    [h,p] = circular_piston(theta_sweep(i),0,fs,fl);
    theta_h(i,:) = h;
end

% radiation pattern
iterPolarplot(freq_sweep, theta_sweep, theta_h, 2, H_size);