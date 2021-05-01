function [AM, AW] = addHead(AW, x, y, z, azimuth, elevation)
% AM: addMic() array
% AW: addWall() array
% x,y,z: coordinates of head in 3D space
%
% Azimuth: angle of receiving sound waves
%          range = [-90;90]º
%
% Elevation: how much the plane xy goes in altitude.
%            range = [-90;90]º
%

% LOAD: head+ears dictionary
%   \_ head wall reflection coefficient = 0
%      (sound normally travels throught the head);
%      Assuming the the head has a diameter of aprox 20cm.
load headnorm mic wall

azimuth = -azimuth;
elevation = -elevation;

% MICROPHONES
for i=1:2
%     % convert ears (<=> microphones) cartesian coordinates
%     % to spherical coordinates
%     [m_th, m_phi, m_r] = cart2sph(mic(i,1), mic(i,2), mic(i,3));
%     % add the function azimuth (<=>theta) (degree to radians)
%     % negative sign to follow my azimuth reference:
%     % left = negative; right = positive
%     m_th = m_th - deg2rad(azimuth);
%     % add the function elevation (<=>phi) (degree to radians in sind())
%     m_phi = m_phi + deg2rad(elevation);
%     % convert spherical coordinates to cartesian coordinates
%     [cart_x, cart_y, cart_z] = sph2cart(m_th,m_phi,m_r);
%     %    \_ updating the cartesians coordinates that
%     %       contain the azimuth and elevation desired
%     mic(i,1:3) = [cart_x+x, cart_y+y, cart_z+z];
    cart_x = mic(i,1);
    cart_z = mic(i,3);
    mic(i,1) = cart_z*sind(elevation) + cart_x*cosd(elevation);
    mic(i,3) = cart_z*cosd(elevation) - cart_x*sind(elevation);
    cart_x = mic(i,1);
    cart_y = mic(i,2);
    mic(i,1) = cart_x*cosd(azimuth) - cart_y*sind(azimuth);
    mic(i,2) = cart_x*sind(azimuth) + cart_y*cosd(azimuth);
end

% HEAD
for i=1:16
%     % convert head cartesian coordinates to spherical coordinates
%     [w_th, w_phi, w_r] = cart2sph(wall(i,1), wall(i,2), wall(i,3));
%     % add the function azimuth (<=>theta) (degree to radians)
%     % negative sign to follow my azimuth reference:
%     % left = negative; right = positive
%     w_th = w_th - deg2rad(azimuth);
%     % add the function elevation (<=>phi) (degree to radians)
%     w_phi = w_phi + deg2rad(elevation);
%     % convert spherical coordinates to cartesian coordinates
%     [cart_x, cart_y, cart_z] = sph2cart(w_th,w_phi,w_r);
%     %    \_ updating the cartesians coordinates that
%     %       contain the azimuth and elevation desired
%     wall(i,1:3) = [cart_x+x, cart_y+y, cart_z+z];
    cart_x = wall(i,1);
    cart_z = wall(i,3);
    wall(i,1) = cart_z*sind(elevation) + cart_x*cosd(elevation);
    wall(i,3) = cart_z*cosd(elevation) - cart_x*sind(elevation);
    cart_x = wall(i,1);
    cart_y = wall(i,2);
    wall(i,1) = cart_x*cosd(azimuth) - cart_y*sind(azimuth);
    wall(i,2) = cart_x*sind(azimuth) + cart_y*cosd(azimuth);
    
    wall(i,:) = wall(i,:) + [x y z];
end



% IMPLEMENTATION: coordinates to objects
%    \_ ears
micdir = [1 -1 0; 1 1 0];

for i=1:2
    cart_x = micdir(i,1);
    cart_z = micdir(i,3);
    micdir(i,1) = cart_z*sind(elevation) + cart_x*cosd(elevation);
    micdir(i,3) = cart_z*cosd(elevation) - cart_x*sind(elevation);
    cart_x = micdir(i,1);
    cart_y = micdir(i,2);
    micdir(i,1) = cart_x*cosd(azimuth) - cart_y*sind(azimuth);
    micdir(i,2) = cart_x*sind(azimuth) + cart_y*cosd(azimuth);
end

mic(1,:) = mic(1,:) + [x y z];
mic(2,:) = mic(2,:) + [x y z];

AM = addMic([], [mic(1,1) mic(1,2) mic(1,3)], @omni, micdir(1,:));
AM = addMic(AM, [mic(2,1) mic(2,2) mic(2,3)], @omni, micdir(2,:));

%    \_ head
for i=[1 5 9 13]
    AW = addWall(AW, wall(i,:), wall(i+1,:), wall(i+2,:), wall(i+3,:), 0);
end
    
end


