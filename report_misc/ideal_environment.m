fclear;clc;close all;
addpath ../simLocUS
addpath ../functions
addpath ../soundfiles/generic
addpath ../structures

%% Read Audio
[snd, fs] = audioread('mosquito.wav');

%% Room Properties
% room prop
room.size = [5 7 3];
room.coeff = [0 0 0];
room.temp = 20;
C = 20.05*sqrt(273.15+room.temp);

%% Frequency and Reflection Properties
MR = 0;
fl = 100;   % lower frequency bound
fh = 20e3;  % upper frequency bound (<=80 kHz due to @KemoL10_TF)

fx = 0:25e3;
fy = zeros(1,length(fx));
fy(fx <= 20e3) = 1;
figure;
semilogx(fx,fy,'b')
xlabel('Frequency (Hz)')
ylabel('Gain')

%% Receiver Properties
% rec prop
rec.loc = [1 3.5 1.7];
rec.struct = 'recstruct.mat';
rec.mic.dmf = [1 2 4 8];
Ndmf = length(rec.mic.dmf);
rec.mic.dist = 0.1*rec.mic.dmf*2;
rec.dir = [0 0];
rec_ = Receiver(rec.struct);

%% Display Room
base_coord = [0 0 0;room.size(1) 0 0; room.size(1) room.size(2) 0;0 room.size(2) 0; 0 0 0];
roof_coord = base_coord;
roof_coord(:,3) = room.size(3);
wall1_coord = [0 0 0;room.size(1) 0 0; room.size(1) 0 room.size(3);0 0 room.size(3); 0 0 0];
wall2_coord = [0 0 0;0 room.size(2) 0; 0 room.size(2) room.size(3);0 0 room.size(3); 0 0 0];
wall3_coord = [0 room.size(2) 0;room.size(1) room.size(2) 0; room.size(1) room.size(2) room.size(3);0 room.size(2) room.size(3); 0 room.size(2) 0];
wall4_coord = [room.size(1) 0 0;room.size(1) room.size(2) 0; room.size(1) room.size(2) room.size(3);room.size(1) 0 room.size(3); room.size(1) 0 0];

f=figure;
fill3(base_coord(:,1),base_coord(:,2),base_coord(:,3), 'k', 'FaceAlpha', 0.5);
hold on
fill3(roof_coord(:,1),roof_coord(:,2),roof_coord(:,3), 'k', 'FaceAlpha', 0.5);
fill3(wall1_coord(:,1),wall1_coord(:,2),wall1_coord(:,3), 'b', 'FaceAlpha', 0.5);
fill3(wall2_coord(:,1),wall2_coord(:,2),wall2_coord(:,3), 'b', 'FaceAlpha', 0.5);
fill3(wall3_coord(:,1),wall3_coord(:,2),wall3_coord(:,3), 'b', 'FaceAlpha', 0.5);
fill3(wall4_coord(:,1),wall4_coord(:,2),wall4_coord(:,3), 'b', 'FaceAlpha', 0.5);
hold off
xlabel('x');
ylabel('y');
zlabel('z');

% system('pdfcrop room_simple.pdf test.pdf');

%% Create Receiver
AW=[];
AM=[];
rec_ = Receiver(rec.struct);
[AM,AW] = addReceiver(AW, rec_, rec.loc, rec.dir);
makeFile('headwall', AW, AM, MR, fs, fl, fh);
displayRoom('headwall');

%% Create Receiver0
AW=[];
AM=[];
rec_ = Receiver(rec.struct);
AM = addReceiver0(rec_, rec.loc, rec.dir);
system('rm headwall.mat');
makeFile('headwall', AW, AM, MR, fs, fl, fh);
displayRoom('headwall');
%%
AW = addDivision(room.size, room.coeff);

points=200;
[sc_xyz, sc_az]= speaker_hcircle2D(room.size, rec.loc, 2, points, 0);


figure
plot(sc_xyz(:,1), sc_xyz(:,2), 'k');
hold on

plot([1 1], [rec.loc(2)-rec.mic.dist(4) rec.loc(2)+rec.mic.dist(4)], 'k--+');
plot(rec.loc(1), rec.loc(2)-0.1*rec.mic.dmf(1),'r>', 'MarkerFaceColor','r');
plot(rec.loc(1), rec.loc(2)+0.1*rec.mic.dmf(1),'r>', 'MarkerFaceColor','r');
hold off
xlim([0 5]);
ylim([0 7]);
xlabel('room x');
ylabel('room y');
pbaspect([1 1 1]);


AZ_DIFF = zeros(points, Ndmf);
DELAY_DIFF = zeros(points, Ndmf);

for i=1:Ndmf
   rangemax = round(rec.mic.dist(i)*fs + 1e-3*fs);
   for j=1:points
        [yL, yR] = sim_stereo(  rec.struct,... 
                                room.size,...
                                rec.loc,...
                                rec.dir,...
                                sc_xyz(j,:),...
                                snd, fs,...
                                1,...
                                rec.mic.dmf(i));

        lr_corr = xcorr(yL,yR, rangemax);
        [value, index] = max(abs(lr_corr));
        delay_index = (index-1) - rangemax;
        delay_t = delay_index*1/fs;
        arg = C*delay_t/rec.mic.dist(i);

        if arg > 1
            arg = 1;
        elseif arg < -1
            arg = -1;
        end

        delta_t = rec.mic.dist(i) * cosd(sc_az(j)) * (1/C);

        AZ_DIFF(j, i) = abs(sc_az(j) - acosd(arg));
        DELAY_DIFF(j, i) = abs(delta_t - delay_t);
   end
   
    f=figure;
    sgtitle('Distance = 2m');
    
    plot(sc_az,AZ_DIFF(:,i));
    xlabel('Source location (º)');
    ylabel('Azimuth error (º)');
    xlim([0 180]);
    
    clc
    fprintf("%.0d%%", round(i/Ndmf*100));
    
    destination="/home/jreis/Nextcloud/Dissertação/microphone_spacing_imgs/";
    path = destination + num2str(i) + ".jpg";
%     exportgraphics(f,path,'Resolution',300)
end