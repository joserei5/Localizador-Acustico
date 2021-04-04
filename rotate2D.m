function rmatrix = rotate2D(cmatrix,azimuth)
% cmatrix: cartesian coordinate matrix [x y]
% azimuth: horizontal rotation in degrees

% cmatrix cartesian coordinates
x = cmatrix(:,1);
y = cmatrix(:,2);

% rotate around z-axis
% (azimuth)
azimuth = - azimuth;
cmatrix(:,1) = x*cosd(azimuth) - y*sind(azimuth);
cmatrix(:,2) = x*sind(azimuth) + y*cosd(azimuth);

rmatrix = cmatrix;
end

