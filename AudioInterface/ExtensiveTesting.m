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
if ispc
    % organize them into a cell without any spaces
    N=size(DIR,1);
    Nf=N-2;
    FILES = cell(N-2,1);
    for i=3:N
        NAME_SPLIT = split(DIR(i,:),' ');
        FILES{i-2,1} = NAME_SPLIT{1};
    end
elseif isunix
    FILES = split(DIR);
    FILES(end) = [];
    Nf = length(FILES);
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
AOA.p2.values = zeros(100 , Nf);
AOA.p2.error = zeros(100 , Nf);

%% Predictor Loop
for i=1:Nf
   % read audio data
   [y,~] = audioread(FILES{i});
   % go through every block
   for j=1:100
    CH.L = y(1+(j-1)*BLK_N:j*BLK_N,1);
    CH.R = y(1+(j-1)*BLK_N:j*BLK_N,2);
    [AOA.p2.values(j,i), ~]= detect_az2(CH,CR,C,REC.d);
   end
end

%% Fix recording "spikes" or "holes"
% and fill the error structure
M2 = median(AOA.p2.values);
DEG_THR = 10; % degree error threshold
THR_MIN2 = M2-DEG_THR;
THR_MAX2 = M2+DEG_THR;

for i=1:Nf
    for j=1:100
        if AOA.p2.values(j,i) > THR_MAX2(i) || AOA.p2.values(j,i) < THR_MIN2(i)
            AOA.p2.values(j,i) = M2(i);
        end
    end
end

% AOA mean
AOA.p2.values = mean(AOA.p2.values);
% AOA error
AOA.p2.error = AOA.ref' - AOA.p2.values;
AOA.p2.error = abs(AOA.p2.error);

%% Algorithm reference-values
% theoretical AOA
r_theta = linspace(0,180,Nf);

% extremes value in seconds
d_t_ext = [(cosd(0)*REC.d)/C (cosd(180)*REC.d)/C];
d_t = linspace(d_t_ext(1), d_t_ext(2), Nf);

% calc AOA with algorithm
arg = C*(d_t+0)./REC.d;
arg(arg<-1) = -1;
arg(arg>1) = 1;
alg_theta = acosd(arg);

%% EXPERIMENTAL Figures
% AOA
mainf=figure;
hold on;
plot(AOA.ref, AOA.p2.values)
plot(AOA.ref, AOA.ref, 'r:')
plot(r_theta, alg_theta, 'g')
hold off;

% AOA error
errorf=figure;
subplot(2,1,1)
plot(AOA.ref, AOA.p2.error)
subplot(2,1,2)
semilogy(AOA.ref, AOA.p2.error)

%% SIMULATOR Figures
% load data
SIM = load('../report_misc/chapter4/sim_line.mat');
% set the minimum error equal to the experimental value
% (pretty log plot)
SIM.AOA_error(SIM.AOA_error==0) = min(AOA.p2.error);

figure(mainf);
hold on;
plot(SIM.TRAJ.th, SIM.AOA,'k')
hold off;

figure(errorf);
subplot(2,1,1)
hold on;
plot(SIM.TRAJ.th, SIM.AOA_error, 'k')
hold off;
subplot(2,1,2)
hold on;
semilogy(SIM.TRAJ.th, SIM.AOA_error, 'k')
hold off;

%% Titles and labels
figure(mainf);
xlabel('Azimuth reference (degrees)')
ylabel('Predicted azimuth (degrees)')
legend( 'Interp (quadr.)',...
        'Azimuth reference', ...
        'Algorithm reference',...
        'Simulator reference');

figure(errorf);
axh=axes(errorf,'visible','off');
axh.Title.Visible='on';
axh.XLabel.Visible='on';
axh.YLabel.Visible='on';
xlabel(axh,'Azimuth reference (degrees)')
ylabel(axh,'Predicted azimuth error (degrees)')