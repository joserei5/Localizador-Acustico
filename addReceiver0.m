function [AM] = addReceiver0(rec, location, angle)
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

% add object
for i=1:size(rec.mic.pos,1)
    if i == 1
        AM = addMic([], rec.mic.pos(i,:), rec.mic.tf, rec.mic.dir(i,:));
    else
        AM = addMic(AM, rec.mic.pos(i,:), rec.mic.tf, rec.mic.dir(i,:));
    end
end