function AOA = getTrajAOA(SRC,ROOM)
% Estimate the AOA from a custom trajectory.
% 
% Input:
%  SRC -> source structure
%  SRC.traj.x: x coordinates for the trajectory points
%  SRC.traj.y: y coordinates for the trajectory points
%
%  REC -> receiver structure
%  ROOM.rec.x:       x coordinate of the center of the receiver
%  ROOM.rec.y:       y coordinate of the center of the receiver
%  ROOM.rec.dx:      microphone spacing
%  ROOM.rec.azimuth: receiver structure azimuth
%
% Output:
%  AOA.theoretical: true values
%  AOA.algorithm:   predicted values
%
% WARNING: it uses the 'recstruct.mat' structure;

% 1) Searching the main folder and getting
% the correct path for the structure file;
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
    error('''Localizador-Acustico/'' not found.');
end
spath = join([repmat('../',1,id),'structures/recstruct.mat']);

% 2) Creating/loading the structure
% Receiver(path,dmf),
% where dmf = distance multiplying factor:
%  |_Usage: 2/diameter * new_diameter/2
%           <=> new_diameter/diameter
%  (1) [0 1 0]:
%      [0 diameter/2 0] * 2/diameter
%  (2) [0 new_diameter/2 0]:
%      [0 1 0] * new_diameter/2
rec_ = Receiver(spath, ROOM.rec.dx/20e-2);
% Update the structure with the correct azimuth
rec_.mic.pos = rotationMatrix(rec_.mic.pos,-ROOM.rec.azimuth,0);

% 3) L and R microphone coordinates
M.L.x = ROOM.rec.x + rec_.mic.pos(1,1);
M.L.y = ROOM.rec.y + rec_.mic.pos(1,2);
M.R.x = ROOM.rec.x + rec_.mic.pos(2,1);
M.R.y = ROOM.rec.y + rec_.mic.pos(2,2);

% 4) Calculate distances:
% 4.a) Point in trajectory to center of receiver
d = sqrt( (ROOM.rec.x - SRC.traj.x).^2 + (ROOM.rec.y - SRC.traj.y).^2 );
% 4.b) Point in trajectory to LEFT microphone
d1 = sqrt( (M.L.x - SRC.traj.x).^2 + (M.L.y - SRC.traj.y).^2 );
% 4.c) Point in trajectory to RIGHT microphone
d2 = sqrt( (M.R.x - SRC.traj.x).^2 + (M.R.y - SRC.traj.y).^2 );

% Calculate theoretical AOA,
% based on a SSS triangle ("side, side, side")
% and applying the law of cosine triangles
%    C
%    /\
%   /  \
%  /    \
%A/______\B
%
% where a = length BC
%       b = length AC
%       c = length AB
%
% so, cos(A) = (b^2 + c^2 − a^2) / (2bc)
%
AOA.theoretical = acosd((ROOM.rec.dx.^2 + d1.^2 - d2.^2)./(2*ROOM.rec.dx*d1));

% Predictor AOA
dd = d1-d2;
arg = dd./ROOM.rec.dx;
arg(arg<-1) = -1;
arg(arg>1) = 1;
AOA.algorithm = acosd(arg);

end

