clear;clc;close all;

addpath ../../simLocUS
addpath ../../functions
addpath ../../soundfiles/generic
addpath ../../structures

debug = true;

%% Igual ao simulador anterior
AW = addWall([], [ 0  0  0], [20  0 0], [20 20 0], [ 0 20 0], 0.5); % chao


AM = addMic([], [-0.1 0 1]);
AM = addMic(AM, [ 0.1 0 1]);

MR = 0;
fs = 48e3;
fl = 0;
fh = 18e3;

makeFile('test_mov', AW, AM, MR, fs, fl, fh);

R = Room();
R.T = 25;
R.H = 30;
R.P = 1.01;

%% Gerar as IRs

% trajectoria linear
v = 10; % velocidade m/s
p_xi = -10; % y inicial
p_xf = 10; % y final
df = v/fs; % salto de trajetoria

p_x = p_xi:df:p_xf;

N = length(p_x);

%%%%%%% Sinal da fonte %%%%%%%%%%%%%%%%
% O numero de amostras de s tem de ser superior a N + a duraçao da resposta
% impulsional

n = 0:N-1;
t = n/fs;

s = 0.1*sin(2*pi*1e3*t) + 0.1*sin(2*pi*4e3*t) + 0.1*sin(2*pi*8e3*t); % sinal emitido pela fonte
s = s'; % tem de ser um vetor coluna

y = zeros(2,N); % 2 microfones e o sinal da saida tera o mesmo numero de pontos que a trajectoria

S = addSpk( [p_x(1) 1 1]); % posiçao inicial da fonte

if debug  
    figure
    pl1 = plot([AM(1).C(1) S.C(1)],[AM(1).C(2) S.C(2)]);
    hold on
    pl2 = plot([AM(2).C(1) S.C(1)],[AM(2).C(2) S.C(2)]);
    plot(AM(1).C(1),AM(1).C(2),'k*')
    plot(AM(2).C(1),AM(2).C(2),'k*')
    ps = plot(S.C(1),S.C(2),'ko');
    axis([p_xi p_xf 0 1.1])
    hold off
    title('Esquema')

    pl1.XDataSource = '[AM(1).C(1) S.C(1)]';
    pl2.XDataSource = '[AM(2).C(1) S.C(1)]';
     ps.XDataSource = 'S.C(1)';
end

tm = 0;
tic
for n=1:N
    
    S.C(1) = p_x(n); % altera a coordenada x da fonte
    
    I = impR('test_mov', S, R); % obtem as novas respostas impulsionais   
    
    M = size(I,2); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% alterado %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n+M-1>N % caso n seja necessario amostras fora do sinal y (assume-se que serao zero) e removem-se da resposta impulsional
        M = N+1-n;
        I(:,M+1:end) = [];
    end
    
    y(:,n:n+M-1) = y(:,n:n+M-1)+s(n).*I;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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