clear;clc;close all;
addpath ../../simLocUS
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../structures

%% Parameters
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
% display auxiliary figures
debug=  true;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
% sound to be reconstructed
wavsample = 'mosquito.wav';
% does it need to be resampled?
doResample = false;
% (if resampling=true) desired sampling frequency
fs2 = 48e3;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
% division size
dSize=  [5 7 3];
% division reflection coefficients
dCoef=  [.7 .8 .5];
% division temperature
rtemp=  20;
% divison humidity
rhum=   55;
% division pressure
rpress= 1.01;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
% order of reflections
MR= 0;
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
% receiver location
recLoc= [1 3.5 1.5];
% receiver direction [azimuth elevation]
recDir= [0 0];
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%% Audio Setup
% load mosquito sample
[s, fs]=    audioread(wavsample);
Ns=         length(s);
t=          (1:Ns)*1/fs;

% resample
if doResample
    s=  resample(s,fs2,fs);
    fs= 48e3;
    Ns= length(s);
    t=  (1:Ns)*1/fs;
end

% % clip audio to 2sec
% sclip = 2/(1/fs);
% s=      s(1:sclip);
% t=      t(1:sclip);

% frequency limits
fl= 100;   % lower frequency bound
fh= 20e3;  % upper frequency bound (<=80 kHz due to @KemoL10_TF)

% plot audio
if debug
    figure
    plot(t,s)
    xlabel('seconds')
    ylabel('Amplitude')
end
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%% Trajectory Setup
% velocity
v=      1;

% first and final point position (y-axis)
p_x=    3;
p_yi=   0.1;
p_yf=   dSize(2)-0.1;

% trajectory jump
df=     v/fs;

% rebuild trajectory
p_y=    p_yi:df:p_yf;
N=      length(p_y);

% initialize stereo channels
y=  zeros(2,N);
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%% Room Setup
% generate division
AW=     addDivision(dSize,dCoef);

% generate head/receiver
rec=    Receiver('recstruct');
AM=     addReceiver0(rec, recLoc, recDir);

% speaker/source initial location
S=  addSpk( [p_x p_y(1) 1]);

% room properties
R=      Room();
R.T=    rtemp;
R.H=    rhum;
R.P=    rpress;

% MAKE divison object
makeFile('test_mov', AW, AM, MR, fs, fl, fh);
%▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%% Simulation
if debug  
    f = figure;
    f.ToolBar = 'none';
    f.MenuBar = 'none';
    f.Resize = 'off';
    pl1 = plot([AM(1).C(1) S.C(1)],[AM(1).C(2) S.C(2)]);
    hold on
    pl2 = plot([AM(2).C(1) S.C(1)],[AM(2).C(2) S.C(2)]);
    plot(AM(1).C(1),AM(1).C(2),'k*')
    plot(AM(2).C(1),AM(2).C(2),'k*')
    ps = plot(S.C(1),S.C(2),'ko');
    pbaspect([1 1 1])
    xlim([0 5]);ylim([0 7]);
    hold off
    title('Esquema')

    pl1.YDataSource = '[AM(1).C(2) S.C(2)]';
    pl2.YDataSource = '[AM(2).C(2) S.C(2)]';
    ps.YDataSource = 'S.C(2)';
end

tm = 0;
tic
for n=1:N
    
    S.C(2) = p_y(n); % altera a coordenada y da fonte
    
    I = impR('test_mov', S, R); % obtem as novas respostas impulsionais   
    
    M = size(I,2); 
    
    if n<M % caso n seja menor que M seria necessário amostras para tempos negativos (assume-se que serão zero) e removem-se da resposta impulsional
        I(:,n+1:end) = [];
        M = n;
    end
    
    y(:,n) = I*s(n:-1:n-M+1); % sinal de saída microfones
    
    if debug   
        ttt = toc;
        if ttt > tm + 1
            refreshdata
            drawnow
            tm = ttt;
            tn = round((N-n)*ttt/(60*n));
            clc
            fprintf('Faltam: %d minutos!\n',tn)
        end
    end
end

spectrogram(y(1,:),2^11,[],[],fs,'yaxis')
ylim([0 10])

soundsc(y,fs)
% audiowrite('mosquito_line.wav',y',fs)