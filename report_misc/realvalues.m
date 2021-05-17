%% Initial cmds
clear;clc;close all;
addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generic
addpath ../structures

%% Data
delta_x = 0.1*2;
room.temp = 20;
C = 20.05*sqrt(273.15+room.temp);

% theoretical AOA
r_theta = 0:0.5:180;

% extremes value in seconds
d_t_ext = [(cosd(0)*delta_x)/C (cosd(180)*delta_x)/C];
d_t = linspace(d_t_ext(1), d_t_ext(2), length(r_theta));

% calc AOA with algorithm
arg = C*(d_t+0)./delta_x;
arg(arg<-1) = -1;
arg(arg>1) = 1;
alg_theta = acosd(arg);
error_theta = r_theta - alg_theta;

%% Plot
% results
figure;
plot(r_theta, r_theta, 'k');
hold on
plot(r_theta, alg_theta, 'r');
hold off
legend('Real AOA', 'Algorithm AOA')
xlabel('source azimuth (º)')
ylabel('output azimuth (º)')

% errors
figure;
plot(r_theta, error_theta, 'b')
xlabel('source azimuth (º)')
ylabel('output azimuth error (º)')