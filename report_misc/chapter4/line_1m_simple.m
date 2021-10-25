clear all;clc;close all;
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../simLocUS
addpath ../../structures

%% Variables
% Room
ROOM.xyz = [3 3 3]; % room coordinates
ROOM.T = 24; % room temperature
% Receiver
REC.xyz = [0.7 ROOM.xyz(2)/2 1.7]; % receiver coordinates
REC.d = 1; % receiver distance in relation to the source
REC.DX = 29.2*1e-2; % receiver: microphone interdistance
REC.th = [-45 45]; % receiver theta
% Algorithm settings
N=9;
MED_t = 1e-3; % maximum estimated delay
C = 20.05*sqrt(273.15+ROOM.T); % sound velocity
AOA = zeros(N,length(REC.th)); % angle of arrival
AOA_error = AOA; % angle of arrival error
TRAJ.th = zeros(N,length(REC.th)); % trajectory theta/azimuth
% Figures
mainf=figure;
errorf=figure;


for i=1:length(REC.th)
%% Prepare theoretical data
    TRAJ.th(:,i) = [45 60:10:120 135]-REC.th(i);
    dm = sqrt( tand(90-TRAJ.th(:,i)).^2 + REC.d );
    d1 = sqrt((0.5*REC.DX - tand(90-TRAJ.th(:,i))).^2 + REC.d);
    d2 = sqrt((0.5*REC.DX + tand(90-TRAJ.th(:,i))).^2 + REC.d);
    dd = d2-d1;
    dt = dd/C;
    
    arg_=               C*dt/(REC.DX);
    if sum(isnan(arg_))>0
        if sum(arg_>0)>0
            arg_(isnan(arg_)) = 1;
        else
            arg_(isnan(arg_)) = -1;
        end
    end
    arg_(arg_>1)=       1;
    arg_(arg_<-1)=      -1;
    AOA(:,i)=           acosd(arg_);
    AOA_error(:,i) = TRAJ.th(:,i) - AOA(:,i);
    AOA_error(:,i) = abs(AOA_error(:,i));

%% Prepare figures
% Create shed and room
    figure(mainf);
    subplot(2,length(REC.th),i)
    % empty room
    r1=rectangle('Position',[0 0 ROOM.xyz(1:2)]);
    % shed
    r2=rectangle('Position',[0 REC.xyz(2)-1/2 1 1]);
    r2.FaceColor = [.9 .9 .9];
    % receiver mark
    hold on;
    plot(REC.xyz(1), REC.xyz(2), 'r*')
    % receiver auxiliary line
    rec_line = [0 -REC.DX/2 0; 0 REC.DX/2 0];
    rec_line = rotationMatrix(rec_line, -REC.th(i), 0);
    rec_line = rec_line + REC.xyz;
    plot(rec_line(:,1), rec_line(:,2), 'r');
    plot(rec_line(:,1), rec_line(:,2), 'r>');
    hold off;
    % adjust axis
    axis equal
    xlim([0 ROOM.xyz(1)])
    ylim([0 ROOM.xyz(2)])
end
AOA_error(AOA_error==0) = 1e-9;


%% FINAL FIGURES
% AOA
    figure(mainf);
    subplot(2,length(REC.th),length(REC.th)+1:2*length(REC.th));
    plot(TRAJ.th, AOA,'b')
    hold on;
    plot(TRAJ.th, TRAJ.th, 'r:')
    hold off;

% AOA error
    figure(errorf);
    subplot(2,1,1)
    plot(TRAJ.th, AOA_error,'b');
    subplot(2,1,2)
    semilogy(TRAJ.th, AOA_error,'b');