clear force all;
clc;
close force all;
addpath ../functions
addpath ../soundfiles/capture/extensive/1m/

%% Audio device info
Fs = 48e3;   % sampling frequency

%% Algorithm Variables
MED_t = 1e-3;                   % maximum estimated delay = 1 ms
MED_N = MED_t * Fs;             % maximum estimated delay (in samples)

ROOM.T = 24.5;                  % room temperature (ºC)
C = 20.05*sqrt(273.15 + ROOM.T);% sound velocity (m/s)
REC.d = 29.2e-2;                % receiver distance = 29.5 cm

CH.fs = Fs;
CR = round(MED_N);     % correlation range=receiver distance + MED (in samples)

%% Prepare variables
% get files from directory
DIR = ls('../soundfiles/capture/extensive/1m/');
% organize them into a cell without any spaces
N=size(DIR,1);
Nf=N-2;
FILES = cell(N-2,1);
for i=3:N
    NAME_SPLIT = split(DIR(i,:),' ');
    FILES{i-2,1} = NAME_SPLIT{1};
end
% sort files numerically
R = cell2mat(regexp(FILES ,'(?<Name>\D+)(?<Nums>\d+)','names'));
tmp = sortrows([{R.Name}' num2cell(cellfun(@(x)str2double(x),{R.Nums}'))]);
FILES = strcat(tmp(:,1) ,cellfun(@(x) num2str(x), tmp(:,2),'unif',0));
for i=1:Nf
    FILES{i} = join([FILES{i},'.wav']);
end
% retrieve theta values from file names
AOA.ref = zeros(Nf,1);
for i=1:Nf
    AOA.ref(i,1)=str2double(regexp(FILES{i},'[\d.]+','match'));
end
% prepare block sizes
BLK_t = 100e-3;             % Block size = 100 ms
BLK_N = BLK_t * Fs;         % 100 blocks of 100ms (10 s of audio)

% prepare the rest of the AOA structure
AOA.p1.values = zeros(100 , Nf);
AOA.p1.error = zeros(100 , Nf);
AOA.p2.values = zeros(100 , Nf);
AOA.p2.error = zeros(100 , Nf);

%% Statistics Loop
for i=1:Nf
   % read audio data
   [y,~] = audioread(FILES{i});
   % go through every block
   for j=1:100
    CH.L = y(1+(j-1)*BLK_N:j*BLK_N,1);
    CH.R = y(1+(j-1)*BLK_N:j*BLK_N,2);
    AOA.p1.values(j,i) = detect_az1(CH,CR,C,REC.d);
    [AOA.p2.values(j,i), ~]= detect_az2(CH,CR,C,REC.d);
   end
end

%% Fix recording "spikes" or "holes"
% and fill the error sstructure
M1 = median(AOA.p1.values);
M2 = median(AOA.p2.values);
DEG_THR = 10; % degree error threshold
THR_MIN1 = M1-DEG_THR;
THR_MAX1 = M1+DEG_THR;
THR_MIN2 = M2-DEG_THR;
THR_MAX2 = M2+DEG_THR;

for i=1:Nf
    for j=1:100
        if AOA.p1.values(j,i) > THR_MAX1(i) || AOA.p1.values(j,i) < THR_MIN1(i)
            AOA.p1.values(j,i) = M1(i);
        end
        if AOA.p2.values(j,i) > THR_MAX2(i) || AOA.p2.values(j,i) < THR_MIN2(i)
            AOA.p2.values(j,i) = M2(i);
        end
        AOA.p1.error(j,i) = AOA.ref(i) - AOA.p1.values(j,i);
        AOA.p2.error(j,i) = AOA.ref(i) - AOA.p2.values(j,i);
    end
end


%% Figures
figure;
plot(AOA.ref, mean(AOA.p1.values))
hold on;
plot(AOA.ref, mean(AOA.p2.values))
plot(AOA.ref, AOA.ref, 'r--')
hold off;
xlabel('Azimuth reference (degrees)')
ylabel('Azimuth predicted (degrees)')
title('Predicted azimuth values')
legend('Not interp.', 'Interp', 'Azimuth reference')

figure;
plot(AOA.ref, abs(mean(AOA.p1.error)))
hold on;
plot(AOA.ref, abs(mean(AOA.p2.error)))
hold off;
xlabel('Azimuth reference (degrees)')
ylabel('Azimuth predicted error (degrees)')
title('Error of predicted azimuth values')
legend('Not interp.', 'Interp')