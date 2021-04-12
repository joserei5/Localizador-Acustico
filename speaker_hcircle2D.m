function [xyz, azimuth] = speaker_hcircle2D(rcoord,hcoord,dist,points,offset)

%                       (!)
% This function follows the implementation of the function
% addDivision(). This works on positive coordinates.

% Room half-circle sound propagation.
% This will create various coordinates for addSpk() function,
% simulating a arc around a semi-cricle (0 to ±180º)
%
%  0º <----------------------(ºvº)----------------------> 180º 
%                 ,agd"Yb
%              ,gdP"
%
%                       _,,dd
%                     ,dP"'
%
%
%                           ,gPP
%                               
%                             \ /
%                              V

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

range1 = 90 + offset;
range2 = -90 + offset;
azimuth = linspace(range1,range2,points);
xyz = zeros(points,3);

xyz(:,1) = H.x+dist*(cosd(azimuth));
xyz(:,2) = H.y-dist*(sind(azimuth));
xyz(:,3) = H.z;

end

