clear force all;
clc;
close force all;

%% Retrieve sample paths
addpath '../../soundfiles/mosquitoes/06_stigma_female'
samples = ls('../../soundfiles/mosquitoes/06_stigma_female');
% windows: remove folder directories and
% .mat file from indexation
if ispc
    samples = samples(3:end-1,:);
else
    return;
end

%% Indexation
Nfiles = size(samples,1);
Nsample = audioinfo(samples(1,:)).TotalSamples;
Ntrack = Nsample*Nfiles; 
% WTrack = zeros(Ntrack,1);
WTrack = [];

%% Compile every sample into a single track
% tic
% for i=1:Nfiles
%     i1 = Nsample*(i-1)+1;
%     i2 = Nsample*i;
%     [y, ~] = audioread(samples(i,:));
%     WTrack(i1:i2,1) = y;
% end
% toc

TMR=tic;
for i=1:Nfiles
    i1 = Nsample*(i-1)+1;
    i2 = Nsample*i;
    [y, ~] = audioread(samples(i,:));
    y = nonzeros(y);
    WTrack = [ WTrack ; y ];
end
toc(TMR);
clear TMR;


%% Figures
fs = audioinfo(samples(1,:)).SampleRate;
% t = (1:Ntrack)*1/fs;
% plot(t,WTrack)
t = (1:length(WTrack))*1/fs;
figure;
plot(t,WTrack)