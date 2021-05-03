function [xyz, azimuth] = speaker_circle2D(rcoord,hcoord,dist,points,offset)

% OUTPUT = [xyz, azimuth]
% 
% [1] rcoord
% room coordinates in WxLxH (input array = [W L H])
%
% [2] hcoord
% head/receiver location coordinates (input array = [x y z])
% 
% [3] dist
% radius of circle (input scalar)
% 
% [4] points
% number of points to represent the circle (input scalar)
% 
% [5] offset
% starting point angle offset (input scalar must be an angle of 0º to 360º)
%
% -------------------------------------------------------------------------
% This function follows the implementation of the function
% addDivision(). This works on positive coordinates.

% Room circle sound propagation.
% This will create various coordinates for addSpk() function,
% simulating a arc around a semi-cricle (0 to ±180º)


% Room struct
R.x = rcoord(1);
R.y = rcoord(2);
R.z = rcoord(3);
% Head struct / Receiver struct
H.x = hcoord(1);
H.y = hcoord(2);
H.z = hcoord(3);

% Adjusting points por even value
if mod(points,2) ~= 0
    points = points + 1;
end
% Clip number of points
if points > 1000
    points = 1000;
elseif points < 10
    points = 10;
end

% Check boundaries
if      (H.x + dist) > R.x
   error('Speaker is outside of boundaries. (x)');
elseif  (H.y + dist) > R.y
   error('Speaker is outside of boundaries. (y)');
end

range1 = 360 + offset;
range2 = 0 + offset;
azimuth = linspace(range1,range2,points);
xyz = zeros(points,3);

xyz(:,1) = H.x+dist*(cosd(azimuth));
xyz(:,2) = H.y-dist*(sind(azimuth));
xyz(:,3) = H.z;

azimuth = flip(azimuth - offset);

end

