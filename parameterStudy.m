clear;clc;close all;
addpath('wav')

room.size = [16 32 3];
room.rec.loc = [8 16 1.7];
room.spk.loc = [4 4 1.7];
room.temp = 20;

rec.struct = 'recstruct.mat';
rec.mic.dmf = [1 2 4 8];
rec.mic.dist = 0.1*rec.mic.dmf*2;
rec.yaw = 0;
rec.pitch = 0;

C = 20.05*sqrt(273.15+room.temp);

points = 500;
distance = linspace(2,8,4);
offset = 90;

% 48000Hz fs
[sndfile, fs] = audioread('cportugal.wav');


%%
D = length(distance);
N = length(rec.mic.dmf);
P = points/2;

EMn = zeros(D,N);

for k=1:D
    [rp, th_azimuth] = speaker_circle2D(   room.size, room.rec.loc,...
                                           distance(k), points, offset);
    th_azimuth = -th_azimuth + offset + 360;

    R1 = zeros(P,N);    % results 1: time delays
    R2 = zeros(P,N);    % results 2: azimuth
    E1 = zeros(P,N);    % error   1: epsilon a.k.a. time delay error

    fprintf("[Distance #%d]", k)
    for i=1:N
        rangemax = round(rec.mic.dist(i)*fs + 1e-3*fs);
        for j=1:P
            [yL, yR] = sim_stereo(  rec.struct,... 
                                    room.size,...
                                    room.rec.loc,...
                                    [rec.yaw, rec.pitch],...
                                    rp(j,:),...
                                    sndfile, fs,...
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

            delta_t = rec.mic.dist(i) * cosd(th_azimuth(j)) * (1/C);
            
            R1(j, i) = delay_t;
            R2(j, i) = acosd(arg);
            E1(j, i) = delay_t - delta_t;
        end
        fprintf("%d ", round(i/N*100))
    end
    disp(" : Completed;");
    
    % plot results
    figure
    subplot(321)
    tstr = "distance of  " + num2str(distance(k)) + "m";
    sgtitle(tstr)
    
    plot(th_azimuth(1:P), R2(:,1), 'r')
    sbtitle = "mic distance of " + num2str(rec.mic.dist(1)) + "m";
    title(sbtitle)
    ylabel('predicted azimuth/yaw (º)')
    xlabel('azimuth/yaw (º)')
    hold on;plot(th_azimuth(1:P), th_azimuth(1:P), 'k--');hold off
    
    subplot(322)
    plot(th_azimuth(1:P), R2(:,2), 'r')
    sbtitle = "mic distance of " + num2str(rec.mic.dist(2)) + "m";
    title(sbtitle)
    ylabel('predicted azimuth/yaw (º)')
    xlabel('azimuth/yaw (º)')
    hold on;plot(th_azimuth(1:P), th_azimuth(1:P), 'k--');hold off
    
    subplot(323)
    plot(th_azimuth(1:P), R2(:,3), 'r')
    sbtitle = "mic distance of " + num2str(rec.mic.dist(3)) + "m";
    title(sbtitle)
    ylabel('predicted azimuth/yaw (º)')
    xlabel('azimuth/yaw (º)')
    hold on;plot(th_azimuth(1:P), th_azimuth(1:P), 'k--');hold off
    
    subplot(324)
    plot(th_azimuth(1:P), R2(:,4), 'r')
    sbtitle = "mic distance of " + num2str(rec.mic.dist(4)) + "m";
    title(sbtitle)
    ylabel('predicted azimuth/yaw (º)')
    xlabel('azimuth/yaw (º)')
    hold on;plot(th_azimuth(1:P), th_azimuth(1:P), 'k--');hold off

    subplot(3,2,[5 6])
    plot(th_azimuth(1:P), E1*1e6)
    sbtitle = "Error on time delay estimation (uS)";
    title(tstr)
    ylabel('epsilon (uS)')
    xlabel('azimuth/yaw (º)')
    legend("mic.dist: " + string(rec.mic.dist))
    
    % estimation errors and statistics
    EMn(k,:) = mean(abs(E1))*1e6;
    
    figure
    tstr = "distance of  " + num2str(distance(k)) + "m";
    sgtitle(["Time resolution  Vs.  Azimuth resolution",tstr])
    for l=1:N
        ang_res = nonzeros(diff(R2(:,l)));
        az_array = transpose(th_azimuth(1:P-1));
        ang_res_mask = and(diff(R2(:,l)), az_array);
        ang_res_az = nonzeros(az_array(ang_res_mask));
        subplot(4,2,l*2-1)
        plot(ang_res_az, ang_res)
        sbtitle = "mic distance of " + num2str(rec.mic.dist(l)) + "m";
        title(sbtitle)
        xlabel("Azimuth (º)")
        ylabel("Resolution")
        
        ang_res = nonzeros(diff(R1(:,l)));
        az_array = transpose(th_azimuth(1:P-1));
        ang_res_mask = and(diff(R1(:,l)), az_array);
        ang_res_az = nonzeros(az_array(ang_res_mask));
%         ang_res = ang_res(2:end-1);
%         ang_res_az = ang_res_az(2:end-1);
        subplot(4,2,l*2)
        plot(ang_res_az, ang_res)
        sbtitle = "mic distance of " + num2str(rec.mic.dist(l)) + "m";
        title(sbtitle)
        xlabel("Azimuth (º)")
        ylabel("Resolution")
%         ytickformat("%.2f")
    end  
end

%%
fprintf("----------------------------------------\n")
fprintf("| d=%.3f | d=%.3f | d=%.3f | d=%.3f|\n", rec.mic.dist)
fprintf("----------------------------------------\n")
for k=1:D
    fprintf("|%.2d |%.2d |%.2d |%.2d|\n", EMn(k,:))
end
fprintf("----------------------------------------\n")

figure
plot(distance, EMn, '.-')
ylabel("Absolute error (uS)")
xlabel("Source distance (m)")