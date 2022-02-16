function ideal_wav(A, T, R, REC)
%ideal_wav.m Creates sound file for the movement of a source in a trajectory
% Inputs:
%   A <=> Audio structure
%   .name   :: file name
%   .fh     :: upper frequency bound
%   .fl     :: lower frequency bound
% 
%   T <=> Trajectory Structure
%   .data
%       .x  :: x points in space
%       .y  :: y points in space
%       .z  :: z points in space
%       .v  :: velocity
%
%   R <=> Room structure
%   .T  :: temperature
%   .H  :: humidity
%   .P  :: pressure
%
%   REC <=> Receiver structure
%   .xyz :: location in space
%   .th  :: theta (amzimuth or yaw)
%   .phi :: phi (elevation or pitch)
%   .DX  :: distance between each microphone
%
% Outputs:
%   local .wav file with completion timestamp
%
% [!] Only works for custom 2D trajectories. Elevation is always constant.


%% Prepare main objects
% floor (no walls)
AW = addWall([], [ 0  0  0], [20  0 0], [20 20 0], [ 0 20 0], 0.5);
% receiver
dir=pwd;
if isunix
    dir=split(dir,'/');
elseif ispc
    dir=split(dir,'\');
end
N=length(dir);
id = -1;
for i=1:N
    if dir{i} == "Localizador-Acustico"
        id = N-i;
        break;
    end
end
if id==-1
    error("Function not located inside 'Localizador-Acustico/functions/'");
end
spath = join([repmat('../',1,id),'structures/recstruct.mat']);
rec_ = Receiver(spath, REC.DX/(2*0.1));
AM = addReceiver0(rec_, REC.xyz, [REC.th REC.phi]);
% load audio
[A.y, A.fs] = audioread(A.name);
% make object
makeFile('app_div_obj', AW, AM, 0, A.fs, A.fl, A.fh);
% room object
R.obj = Room();
R.obj.T = R.T;
R.obj.H = R.H;
R.obj.P = R.P;

%% Prepare trajectory (with interpolation)
% generate waypoints
waypoints = [T.data.x, T.data.y, T.data.z];
% trajectory differences
tt.DX = diff(T.data.x);
tt.DX = [0; tt.DX];
tt.DY = diff(T.data.y);
tt.DY = [0; tt.DY];
% trajectory distance estimation
tt.D = sqrt(tt.DX.^2 + tt.DY.^2);
% trajectory time
tt.T = tt.D ./ T.data.v;
tt.T = cumsum(tt.T);
% define waypointTrajectory parameters
SRATE = A.fs;
FRAMES = SRATE;
% waypointTrajectory setup
trajectory = waypointTrajectory(waypoints,...
                                'TimeOfArrival',tt.T,...
                                'SampleRate', SRATE,...
                                'SamplesPerFrame', FRAMES);
% initialize dumping arrays
tx = []; ty = []; tz = [];
% get last dump index
last_dump = round(round(tt.T(end))*(SRATE/FRAMES));
% initialize time estimation variables
tic;tm=0;tn=0;

% initialize informative message box
mb = msgbox("? minutes remaining","[ideal_wav.m: Trajectory Interpolation");
% start the waypointTrajectory dumping cycle
for i=1:last_dump
   dump = trajectory();
   tx = [tx; dump(:,1)];
   ty = [ty; dump(:,2)];
   tz = [tz; dump(:,3)];

   % tn = minutes remaining
    ttt = toc;
    if ttt > tm + 1
        tm = ttt;
        tn = round((last_dump-i)*ttt/(60*i));
    end
   txt = sprintf("%d minutes remaining",tn);
   set(findobj(mb,'Tag','MessageBox'),'String',txt);
   drawnow;
end
% close message box
delete(mb);

%% Generate IR
% initialize stereo channels
N = length(tx);
y = zeros(2,N);

% check if sound samples do not cover trajectory
% (if they don't: repeat the samples until the length
% of the array is the same length of the interpolated
% trajectory array)
if N > length(A.y)
    A.y = repmat(A.y, ceil(N/length(A.y)), 1);
end

% speaker/source initial location
S = addSpk([tx(1) ty(1) tz(1)]);

% initialize time estimation variables
% time variables
tm = 0; tn=0;
% start timer/clock
clk = tic;

% RECONSTRUCTION LOOP
% initialize informative message box
mb = msgbox("? minutes remaining","[ideal_wav.m: Spatial reconstruction");
for n=1:N
    % skip on: out-of-boundaries coordinates
    if tx(n) <= 0
        continue
    end
    if ty(n) <= 0
        continue
    end
    if tz(n) <= 0
        continue
    end
    % fix NaN numbers
    if n>1
        if isnan(tx(n))
            tx(n) = tx(n-1);
        end
        if isnan(ty(n))
            ty(n) = ty(n-1);
        end
        if isnan(tz(n))
            tz(n) = tz(n-1);
        end
    else
        if isnan(tx(n))
            tx(n) = T.data.x(1);
        end
        if isnan(ty(n))
            ty(n) = T.data.y(1);
        end
        if isnan(tz(n))
            tz(n) = T.data.z(1);
        end
    end

    % feed data to speaker/moving object
    S.C(1) = tx(n);
    S.C(2) = ty(n);
    S.C(3) = tz(n);

    % impulse responses
    I = impR('app_div_obj', S, R.obj);
    M = size(I,2); 

    % when samples outside of signal y
    % are not needed -
    % (assuming that they are null values)
    % - those samples are removed from
    % the impulse response I
    if n+M-1>N
        M = N+1-n;
        I(:,M+1:end) = [];
    end
    y(:,n:n+M-1) = y(:,n:n+M-1)+A.y(n).*I;

    % tn = remaining minutes
    ttt = toc(clk);
    if ttt > tm + 1
        tm = ttt;
        tn = round((N-n)*ttt/(60*n));
    end
    
   txt = sprintf("%d minutes remaining",tn);
   set(findobj(mb,'Tag','MessageBox'),'String',txt); 
   drawnow;
end
% close message box
delete(mb);

% write stereo sample
% (same name as the trajectory)
path = join([repmat('../',1,id),'soundfiles/generated/']);
fname = datestr(now,'ddmmyyyy_HHMMSS');
fpath = join([path,fname,'.wav']);
audiowrite(fpath,y',A.fs)

