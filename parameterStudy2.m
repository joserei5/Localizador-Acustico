%% Initial cmds
clc
clear
close all

%% Parameters
N = 40;
delta_x = transpose(linspace(0.05,20,N));
epsilon = transpose(linspace(-500*1e-6,500*1e-6,N));
% distances
d = transpose(linspace(2,20,N));

theta = 15;

room.temp = 20;
C = 20.05*sqrt(273.15+room.temp);

%% Theoretical SURF plot
% source P(x,y) coordinates
source = [d*cosd(theta) d*sind(theta)];

sensor1 = [-delta_x/2 zeros(N, 1)];
sensor2 = [+delta_x/2 zeros(N, 1)];

% visualize source P(x,y)
figure;
plot(0,0,'b+');hold on
plot(source(:,1), source(:,2), 'k*');hold off

% calculate distance to sensor1 and sensor2
d1 = sqrt((sensor1(:,1)-source(:,1)).^2 + (sensor1(:,2)-source(:,2)).^2);
d2 = sqrt((sensor2(:,1)-source(:,1)).^2 + (sensor2(:,2)-source(:,2)).^2);

% difference of distance in sensor1 in relation to sensor2
delta_d = d1 - d2;
delta_t = delta_d/C;

% prepare variables for 3D curve
[dX, d] = meshgrid(delta_x, d);
% apply algorithm
arg = C*(delta_t+0)./dX;
arg(arg<-1) = -1;
arg(arg>1) = 1;
new_azimuth = acosd(arg);
error_az = ones(N,1)*theta - new_azimuth;

% visualize curve f(delta_x, distance, error)
figure
h = surf(dX, d, error_az, 'FaceColor', 'flat');
% set(h,'LineStyle','none')
title("Epsilon = 0uS")
xlabel('microphone spacing (m)')
ylabel('source distance (m)')
zlabel('azimuth error (º)')