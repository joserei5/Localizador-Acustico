clear;clc;close all;
addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generic
addpath ../structures

%% Audio Properties
[inAudio.s, inAudio.fs] = audioread('mosquito.wav');
inAudio.fl = 100;   % lower frequency bound
inAudio.fh = 20e3;  % upper frequency bound (<=80 kHz due to @KemoL10_TF)

%% Division/Room Properties
% basic parameters
div.size = [5 7 3];
div.coeff = [0 0 0];

% reflection Properties
div.MR = 0;

% generate Room Structure
div.R = Room();
div.R.T = 20;
div.R.H = 30;
div.R.P = 1.01;

% sound velocity
C = 20.05*sqrt(273.15+div.R.T);

%% Receiver Properties
rec.type = 0;
rec.loc = [1 3.5 1.7];
rec.struct = 'recstruct.mat';
rec.mic.dmf = 1;
rec.mic.dist = 0.1*rec.mic.dmf*2;
rec.th = 0;
rec.phi = 0;

%% Speaker Properties
% distance points
d.pts = 100;

% IRREGULAR trajectory along y-axis
traj_x_offset = +3;
spk.loc = zeros(d.pts, 3);
spk.loc(:,3) = rec.loc(3)*ones(d.pts,1);
spk.loc(:,2) = linspace(0.1, div.size(2)-0.1, d.pts);
spk.loc(:,1) = rec.loc(1)*ones(d.pts,1) + traj_x_offset + ...
               transpose(0.5*sin(linspace(0,18*pi,d.pts))) + randn(d.pts,1)/10;

%% Plots
traj_f = figure;

% room
rectangle('Position',[0 0 div.size(1) div.size(2)])
hold on;

% trajectory
plot(spk.loc(:,1), spk.loc(:,2), 'b')
plot(spk.loc(end,1), spk.loc(end,2), 'b^')

% receivers + center
rec_offset = load(rec.struct); rec_offset=rec_offset.rec.mic.pos(1,2);
plot(rec.loc(1), rec.loc(2), 'r+')
plot(rec.loc(1), rec.loc(2)+rec_offset, 'r>')
plot(rec.loc(1), rec.loc(2)-rec_offset, 'r>')
spk_ptx = spk.loc(1,1);
spk_pty = spk.loc(1,2);
traj_pt = plot(spk_ptx, spk_pty, 'k+');
traj_pt.XDataSource = 'spk_ptx';
traj_pt.YDataSource = 'spk_pty';

% detection line
det_src_x = spk.loc(1,1);
det_src_y = spk.loc(1,2);
det_line_x = [rec.loc(1) det_src_x];
det_line_y = [rec.loc(2) det_src_y];
det_plt= plot(det_line_x, det_line_y, 'r');
det_plt.XDataSource = 'det_line_x';
det_plt.YDataSource = 'det_line_y';

% theoretical azimuth line
taz_src_x = spk.loc(1,1);
taz_src_y = spk.loc(1,2);
taz_line_x = [rec.loc(1) taz_src_x];
taz_line_y = [rec.loc(2) taz_src_y];
taz_plt= plot(taz_line_x, taz_line_y, 'k');
taz_plt.XDataSource = 'taz_line_x';
taz_plt.YDataSource = 'taz_line_y';

% display and correct axes
axis equal
xlim([0 div.size(1)])
ylim([0 div.size(2)])
drawnow

%% Generate Stereo Samples
% and apply detection algorithm

% initialize Angle-of-arrival array (azimuth values)
AOA = zeros(d.pts,1);

% determine theoretical value of AOA
AOA_t = 180-atan2d(spk.loc(:,1)-rec.loc(1), spk.loc(:,2)-rec.loc(2));

% update correlation range
CR = round(rec.mic.dist + 1e-3*inAudio.fs);

% initialize channel structure
CH.L = 0;
CH.R = 0;
CH.fs = inAudio.fs;

for i=1:d.pts
    % refresh speaker points
    spk2.loc = spk.loc(i,:);
    
    % create stereo samples
    [CH.L, CH.R] = sim_stereo(rec, div, spk2, inAudio, 1);
    
    % azimuth algorithm
    AOA(i) = detect_az1(CH, CR, C, rec.mic.dist);
    
    % refresh trajectory pointer plot (traj_pt)
    spk_ptx = spk2.loc(1);
    spk_pty = spk2.loc(2);
    
    det_d = sqrt((rec.loc(1)-spk2.loc(1))^2+(rec.loc(2)-spk2.loc(2))^2);
    det_src_x = rec.loc(1)+det_d*cosd(AOA(i)-90);
    det_src_y = rec.loc(2)+det_d*sind(AOA(i)-90);
    det_line_x = [rec.loc(1) det_src_x];
    det_line_y = [rec.loc(2) det_src_y];
    
    taz_src_x = rec.loc(1)+det_d*cosd(AOA_t(i)-90);
    taz_src_y = rec.loc(2)+det_d*sind(AOA_t(i)-90);
    taz_line_x = [rec.loc(1) taz_src_x];
    taz_line_y = [rec.loc(2) taz_src_y];
    
    refreshdata;
    drawnow;
    
    % progress status bar
    if mod(i,10)==0
        fprintf(repmat('\b', 1, 25));
        fprintf("Stereo + Detection: %.0d%%\n", round(i/d.pts*100))
    end
end

% compare results
results_p = figure;
plot(spk.loc(:,2),AOA,'r')
hold on
plot(spk.loc(:,2),AOA_t,'k')
xlabel('Y-axis position (m)')
ylabel('Theta / Azimuth (º)')
legend('algorithm AOA','correct AOA')

% theta resolution
AOA_s = [0; diff(AOA)];
AOA_s_unique = AOA_s(AOA_s~= 0);
AOA_sl = length(AOA_s_unique);

s_v = 1;
s_flag = 0;
for i=1:d.pts
    if AOA_s(i) == 0
        AOA_s(i) = AOA_s_unique(s_v);
        s_flag = 1;
    else
        if s_flag==1
            s_v = s_v+1;
        end
        s_flag = 0;
    end
end

th_step_p = figure;
t = tiledlayout(1,1);
ax1 = axes(t);
ax2 = axes(t);
plot(ax1, spk.loc(:,2),AOA_s,'b')
xlabel(ax1,'Y-axis position (m)')
ylabel(ax1,'Theta resolution (º)')
plot(ax2, AOA_t ,AOA_s,'b')
ax2.XAxisLocation = 'top';
ax2.YTick = [];
ax2.XTick = round(AOA_t(1),-1):10:round(AOA_t(end),-1);
xlabel(ax2,'Azimuth / Theta (º)')


% theta detection error
AOA_e = AOA_t - AOA;
th_error_p = figure;
t = tiledlayout(1,1);
ax1 = axes(t);
ax2 = axes(t);
plot(ax1, spk.loc(:,2),AOA_e,'b')
xlabel(ax1,'Y-axis position (m)')
ylabel(ax1,'Theta error (º)')
plot(ax2, AOA_t ,AOA_e,'b')
ax2.XAxisLocation = 'top';
ax2.YTick = [];
ax2.XTick = round(AOA_t(1),-1):10:round(AOA_t(end),-1);
xlabel(ax2,'Azimuth / Theta (º)')