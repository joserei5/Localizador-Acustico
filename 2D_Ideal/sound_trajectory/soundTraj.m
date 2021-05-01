clear;clc;close all;
addpath ../../simLocUS
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../structures

%% INIT
% import audio
[y, fs] = audioread('mosquito.wav');

% plot mosquito sample
t = (1:length(y))*1/fs;
figure;
plot(t*1e3,y)
xlabel('t (ms)')
ylabel('Amplitude')

% room properties
room.size = [5 7 3];
room.rec.loc = [2.5 3.5 1.7];
room.temp = 20;

% sound velocity properties
C = 20.05*sqrt(273.15+room.temp);

% plot room 2D (receiver)
f = figure;
f.ToolBar = 'none';
f.MenuBar = 'none';
f.Resize = 'off';
f.Position(3) = 400;
f.Position(4) = 400;
rectangle('Position',[0 0 room.size(1) room.size(2)])
hold on;
plot(room.rec.loc(1), room.rec.loc(2),'b+')
hold off;
text(room.rec.loc(1)-0.35, room.rec.loc(2)-0.2,'Receiver', 'Color', 'b')

% receiver properties
rec.struct = 'recstruct.mat';
rec.mic.dmf = 1;
rec.mic.dist = 0.1*rec.mic.dmf*2; %[!] modify struct so i can change diam of receiver [!]
rec.yaw = 0;
rec.pitch = 0;

% define maximum range for xcorrelation
rangemax = round(rec.mic.dist+1e-3*fs);

% half-circle trajectory
distance = 1.5;
points = 1000;
[xyz, azimuth] = speaker_hcircle2D(  room.size, room.rec.loc,    ...
                                    distance, points, 0);
text(xyz(1,1), xyz(1,2)-0.5, 'Start')
                                
% plot room 2D (speakers)
hold on;
plot(xyz(:,1), xyz(:,2), 'k-')
hold off;


%% GEN TRAJ STEREO
% predicted azimuth/yaw array
azimuth = zeros(points,1);

% determine samples length
Ny = length(y);

% clip sound into <points> blocks
samples_c = floor(Ny/points);

% fade in + fade out filter
filt_f = figure;
fth1 = round(samples_c*0.3);
fth2 = round(samples_c*0.15);
fadein = sin(linspace(0,pi/2,fth1));
fadeout = cos(linspace(0,pi/2,fth2));
nofade = linspace(1,1,samples_c-fth1-fth2);
wavfilter = [fadein nofade fadeout];
plot(1:samples_c, wavfilter)
xlim([0 samples_c])

y_clip = zeros(samples_c, points);
for i=1:points
    i_start = 1+samples_c*(i-1);
    i_end   = samples_c*i;
%     y_clip(:,i) = y(i_start:i_end,1).*wavfilter';
    y_clip(:,i) = y(i_start:i_end,1);
end

% determine how many samples are considered to
% simulate 10% of the time of each window
sth = samples_c*(1/fs)*0.1;
s10p = length(t(t<sth));
% determine samples of non-overlapped window
Ncut = samples_c-s10p;
% trajectory sound total samples
N = Ncut*points;

% trajectory sound array with 2 channels (L and R)
traj_snd = zeros(N,2);

for i=1:points
    [L, R] =    sim_stereo( rec.struct, room.size, room.rec.loc,    ...
                [rec.yaw, rec.pitch], xyz(i,:),                     ...
                y_clip(:,i), fs, 1);

    i_start = 1+Ncut*(i-1);
    i_end   = Ncut*i;
    traj_snd(i_start:i_end,1) = L(1:Ncut,1);
    traj_snd(i_start:i_end,2) = R(1:Ncut,1);
        
    lr_corr = xcorr(L,R, rangemax);
    [value, index] = max(abs(lr_corr));
    delay_index = (index-1) - rangemax;
    delay_t = delay_index*1/fs;

    arg = C*delay_t/rec.mic.dist;
    arg(arg<-1) = -1;
    arg(arg>1) = 1;
    azimuth(i,1) = acosd(arg);
    
    fprintf("%.0d%%", round(i/points*100))
    clc
end

sound(traj_snd, fs)
% figure(f)
% hold on;
% for i=1:points
%     plot(xyz(i,1), xyz(i,2), 'r.')
%     drawnow
%     pause(Ncut*1/fs-0.1*(Ncut*1/fs))
% end