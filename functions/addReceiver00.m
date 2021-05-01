function [AM] = addReceiver00(rec, location, angle)
% AM: addMic() array
% AW: addWall() array
%
% location: [x y z] coordinates of head in 3D space
%
% angle = [theta phi]
% theta = azimuth = yaw: range = [-180;180]º
% phi = elevation = pitch: range = [-90;90]º
%
% rec = structure containing microphone & walls/faces
% cartesian coordinates

%% INIT
% make positive yaw, i.e., right side
yaw = -angle(1,1);
% make positive pitch, i.e., upwards
pitch = -angle(1,2);


%% IMPLEMENTATION OF MICROPHONES
% define direction
rec.mic.dir = rotationMatrix(rec.mic.dir,yaw,pitch);

% rotate microphone to desired angles
rec.mic.pos = rotationMatrix(rec.mic.pos,yaw,pitch);

% translate microphone to desired location
rec.mic.pos = rec.mic.pos + location;

% add ONLY 1 object
AM = addMic([], rec.mic.pos(1,:), rec.mic.tf, rec.mic.dir(1,:));