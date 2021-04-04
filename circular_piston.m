function [h,p] = circular_piston(theta, phi, fs, fl, fh)
% Theta and phi -> MATLAB spheric coordinates
% Considering a piston disc radius of ~0.5cm
% a = 5e-3;
a = 0.5*10^(-2);
% Besselj = bessel function of first kind.

% theta = pointer ^r angle in relation to x-axis (see sph2cart)
% pointer vector r^
% r_hat = [cosd(theta)*cosd(phi) sind(phi) sind(theta)*cosd(phi)];
[r_x r_y r_z] = sph2cart(theta,phi,1);
r_hat = [r_x r_y r_z];

% versor vector x
x_hat = [1 0 0];
% theta angle of the circular piston
% p_theta = acos(r_hat * x_hat')/(norm(r_hat)*norm(x_hat));
p_theta = acos(dot(r_hat,x_hat));
% disp(p_theta)

% Create equispaced frequency points until fs/2
n = 1000;
f = linspace(0,fs/2,n);

% Impulse response of the circular piston
lambda = 340./f;
arg = ((2*pi*a)./lambda)*sin(p_theta);
h = 2*besselj(1,arg)./arg;
% figure;
% plot(f,h)
% xlabel('f (Hz)')
% ylabel('J1')
% title(' J1(f) ')

% Exterior case of the circular piston
K = 0.2;
h = h*(1-K+K*cos(p_theta));
%h = [h zeros(1,n)]; % fs/2 to fs with 0's
p = 0;