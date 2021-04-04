function [AW] = addDivision(cmatrix,rcmatrix)

%             +-----------------+
%            /|                 |
%           / |                 |
%          /  |                 |
%         /   +---------------- +
%        /   /                 /
%       +   /                 /
%       |  /                 / Y
%     Z | /                 /
%       |/                 /
%       +-----------------+
%               X
%
% cartesian coordinates [x y z]
x = cmatrix(1,1);
y = cmatrix(1,2);
z = cmatrix(1,3);
%
% rcmatrix = [wall ceil floor] reflection coefficients
wall_rc = rcmatrix(1,1);
ceil_rc = rcmatrix(1,2);
floor_rc = rcmatrix(1,3);
%
% AW = addWall() array
% ceiling
AW = addWall([], [0 0 z], [x 0 z], [x y z], [0 y z], ceil_rc);
% floor
AW = addWall(AW, [0 0 0], [x 0 0], [x y 0], [0 y 0], floor_rc);
% walls
AW = addWall(AW, [0 0 0], [0 0 z], [x 0 z], [x 0 0], wall_rc);
AW = addWall(AW, [0 0 0], [0 0 z], [0 y z], [0 y 0], wall_rc);
AW = addWall(AW, [x 0 0], [x 0 z], [x y z], [x y 0], wall_rc);
AW = addWall(AW, [0 y 0], [0 y z], [x y z], [x y 0], wall_rc);

end

