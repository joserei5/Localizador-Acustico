function D=orthoNVec(Dx)
%
%  orthoNVec - Compute a set of orthonormal vectors.
%
%  D = orthoNVec(Dx)
%
%  Input:
%   Dx - direction of the vector Dx.
%
%  Copyright 2011 LocUS Project.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 2.0$    $Date: 2011/11/18 16:42$

    Dx = Dx./norm(Dx);

    th1 = (180/pi)*atan2(Dx(3),Dx(2));

    v2y = Dx(2)*cosd(-th1)-Dx(3)*sind(-th1);

    th2 = (180/pi)*atan2(v2y,Dx(1));

    Dx = [cosd(th2); cosd(th1)*sind(th2); sind(th1)*sind(th2)];
    Dy = [-sind(th2); cosd(th1)*cosd(th2); sind(th1)*cosd(th2)];
    Dz = [0; -sind(th1); cosd(th1)];
    
    D = [Dx, Dy, Dz];
end