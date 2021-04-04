function rmatrix = rotationMatrix(cmatrix,azimuth,elevation)
% cmatrix: cartesian coordinate matrix [x y z]
% azimuth: horizontal rotation in degrees
% elevation: vertical rotation in degrees

% cmatrix cartesian coordinates
x = cmatrix(:,1);
y = cmatrix(:,2);
z = cmatrix(:,3);

% rotate around y-axis
% (elevation)
cmatrix(:,1) = z*sind(elevation) + x*cosd(elevation);
cmatrix(:,3) = z*cosd(elevation) - x*sind(elevation);

% retrieve updated coordinates
x = cmatrix(:,1);

% rotate around z-axis
% (azimuth)
cmatrix(:,1) = x*cosd(azimuth) - y*sind(azimuth);
cmatrix(:,2) = x*sind(azimuth) + y*cosd(azimuth);

rmatrix = cmatrix;

end

