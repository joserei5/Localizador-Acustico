function [AM, AW] = addHumanHead(AW,x,y,z,direction,azimuth)
% AM: addMic() array
% AW: addWall() array
% x,y,z: coordinates of head in 3D space
% direction: axis where head is looking;
% ('±xx','±yy')
% ('±zzx', '±zzy') ---> pointing in zz axis but with ears built in
%                       a specific xx/yy axis
%                  +zz
%                   | -xx
%                   | / 
%                   |/
%        -yy <------+------> +yy
%                  /|
%                 / |
%               +xx |
%                  -zz
%
% Azimuth: angle of receiving sound waves[-90º;90º]
%                     0º
%                ~         ~
%             ~               ~
%          ~                     ~
%         ~                       ~
%   -90º <----------(o^o)----------> 90º 

r_angle = cosd(azimuth);
l_angle = -cosd(azimuth);

% head wall reflection coefficient = 0
%(sound normally travels throught the head)
% Assuming the the head has a diameter of aprox 20cm

switch direction
    case '+xx'
        AM = addMic([], [x y-0.1 z], @omni, [1 l_angle 0]);
        AM = addMic(AM, [x y+0.1 z],  @omni, [1 r_angle 0]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z-0.1],...
                    [x+0.1 y-0.1 z+0.1], [x-0.1 y-0.1 z+0.1], 0);
        %right head wall 
        AW = addWall(AW, [x-0.1 y+0.1 z-0.1], [x+0.1 y+0.1 z-0.1],...
                    [x+0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z+0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x+0.1 y+0.1 z-0.1], [x+0.1 y+0.1 z+0.1],...
                    [x+0.15 y z+0.1], [x+0.15 y z-0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x+0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z+0.1],...
                    [x+0.15 y z+0.1], [x+0.15 y z-0.1], 0);
    case '-xx'
        AM = addMic([], [x y-0.1 z], @omni, [-1 l_angle 0]);
        AM = addMic(AM, [x y+0.1 z],  @omni, [-1 r_angle 0]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z-0.1],...
                    [x+0.1 y-0.1 z+0.1], [x-0.1 y-0.1 z+0.1], 0);
        %right head wall
        AW = addWall(AW, [x-0.1 y+0.1 z-0.1], [x+0.1 y+0.1 z-0.1],...
                    [x+0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z+0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x-0.1 y+0.1 z-0.1], [x-0.1 y+0.1 z+0.1],...
                    [x-0.15 y z-0.1], [x-0.15 y z+0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x-0.1 y-0.1 z+0.1],...
                    [x-0.15 y z+0.1], [x-0.15 y z-0.1], 0);
    case '+yy'
        AM = addMic([], [x-0.1 y z], @omni, [l_angle 1 0]);
        AM = addMic(AM, [x+0.1 y z],  @omni, [r_angle 1 0]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x-0.1 y-0.1 z+0.1],...
                    [x-0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z-0.1], 0);
        %right head wall
        AW = addWall(AW, [x+0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z+0.1],...
                    [x+0.1 y+0.1 z+0.1], [x+0.1 y+0.1 z-0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x-0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z-0.1],...
                    [x y+0.15 z-0.1], [x y+0.15 z+0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x+0.1 y+0.1 z+0.1], [x+0.1 y+0.1 z-0.1],...
                    [x y+0.15 z-0.1], [x y+0.15 z+0.1], 0);
    case '-yy'
        AM = addMic([], [x-0.1 y z], @omni, [l_angle -1 0]);
        AM = addMic(AM, [x+0.1 y z],  @omni, [r_angle -1 0]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x-0.1 y-0.1 z+0.1],...
                    [x-0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z-0.1], 0);
        %right head wall
        AW = addWall(AW, [x+0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z+0.1],...
                    [x+0.1 y+0.1 z+0.1], [x+0.1 y+0.1 z-0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x-0.1 y-0.1 z+0.1], [x-0.1 y-0.1 z-0.1],...
                    [x y-0.15 z-0.1], [x y-0.15 z+0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x+0.1 y-0.1 z+0.1], [x+0.1 y-0.1 z-0.1],...
                    [x y-0.15 z-0.1], [x y-0.15 z+0.1], 0);
    case '+zzx'
        AM = addMic([], [x y-0.1 z], @omni, [0 l_angle 1]);
        AM = addMic(AM, [x y+0.1 z],  @omni, [0 r_angle 1]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z-0.1],...
                    [x+0.1 y-0.1 z+0.1], [x-0.1 y-0.1 z+0.1], 0);
        %right head wall
        AW = addWall(AW, [x-0.1 y+0.1 z-0.1], [x+0.1 y+0.1 z-0.1],...
                    [x+0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z+0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x+0.1 y z+0.15], [x-0.1 y z+0.15],...
                    [x-0.1 y-0.1 z+0.1], [x+0.1 y-0.1 z+0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x+0.1 y z+0.15], [x-0.1 y z+0.15],...
                    [x-0.1 y+0.1 z+0.1], [x+0.1 y+0.1 z+0.1], 0);
    case '-zzx'
        AM = addMic([], [x y-0.1 z], @omni, [0 l_angle -1]);
        AM = addMic(AM, [x y+0.1 z],  @omni, [0 r_angle -1]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z-0.1],...
                    [x+0.1 y-0.1 z+0.1], [x-0.1 y-0.1 z+0.1], 0);
        %right head wall
        AW = addWall(AW, [x-0.1 y+0.1 z-0.1], [x+0.1 y+0.1 z-0.1],...
                    [x+0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z+0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x+0.1 y z-0.15], [x-0.1 y z-0.15],...
                    [x-0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z-0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x+0.1 y z-0.15], [x-0.1 y z-0.15],...
                    [x-0.1 y+0.1 z-0.1], [x+0.1 y+0.1 z-0.1], 0);
    case '+zzy'
        AM = addMic([], [x-0.1 y z], @omni, [l_angle 0 1]);
        AM = addMic(AM, [x+0.1 y z],  @omni, [r_angle 0 1]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x-0.1 y-0.1 z+0.1],...
                    [x-0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z-0.1], 0);
        %right head wall
        AW = addWall(AW, [x+0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z+0.1],...
                    [x+0.1 y+0.1 z+0.1], [x+0.1 y+0.1 z-0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x-0.1 y+0.1 z+0.1], [x y+0.1 z+0.15],...
                    [x y-0.1 z+0.15], [x-0.1 y-0.1 z+0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x+0.1 y+0.1 z+0.1], [x y+0.1 z+0.15],...
                    [x y-0.1 z+0.15], [x+0.1 y-0.1 z+0.1], 0);
    case '-zzy'
        AM = addMic([], [x-0.1 y z], @omni, [l_angle 0 -1]);
        AM = addMic(AM, [x+0.1 y z],  @omni, [r_angle 0 -1]);
        %left head wall
        AW = addWall(AW, [x-0.1 y-0.1 z-0.1], [x-0.1 y-0.1 z+0.1],...
                    [x-0.1 y+0.1 z+0.1], [x-0.1 y+0.1 z-0.1], 0);
        %right head wall
        AW = addWall(AW, [x+0.1 y-0.1 z-0.1], [x+0.1 y-0.1 z+0.1],...
                    [x+0.1 y+0.1 z+0.1], [x+0.1 y+0.1 z-0.1], 0);
        %left orthogonal wall
        AW = addWall(AW, [x-0.1 y+0.1 z-0.1], [x y+0.1 z-0.15],...
                    [x y-0.1 z-0.15], [x-0.1 y-0.1 z-0.1], 0);
        %right orthogonal wall
        AW = addWall(AW, [x+0.1 y+0.1 z-0.1], [x y+0.1 z-0.15],...
                    [x y-0.1 z-0.15], [x+0.1 y-0.1 z-0.1], 0);
    otherwise
        tip1 = 'Try with ''+xx'',''-xx'',''+yy'', ''-yy'',';
        tip2 = '''+zzx'',''-zzx'',''+zzy'', ''-zzy''.';
        error('%s is not a valid direction.\n%s\n%s',direction,tip1,tip2)
end

end

