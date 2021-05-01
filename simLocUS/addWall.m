function AW=addWall(AW,P1,P2,P3,P4,RCf)
%
%  addWall - Add a surface (Wall) to the surface array AW.
%
%  AW = addWall(AW, P1, P2, P3, P4, tR)
%
%   INPUT:
%   AW    - Surface array that you want to add a new surface;
%           (first call AW = [])
%   Px    - The four vertices of the surface
%           (It assumes that the surface is a quadrilateral);
%   RCf   - Reflection coefficient
%           or
%           Reflection coefficient function handle (reflection impulse response);
%           It must be a function of type [h,p] = func(fs,fL,fH); where h is
%           the impulse response of the wall and p the position of the first
%           impulse response sample; fs the sampling frequency, fL and fH
%           are the lower and upper frequency bounds for band optimization
%
%  Copyright 2011 LocUS Project.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 2.1$    $Date: 2011/12/21 15:50$

if(nargin<4 || nargin>6 || nargin==5)
    error('Invalid input arguments.');
end

if(nargin==4)
    x = AW;
    y = P1;
    z = P2;
    
    AW = struct('P1',{},'P2',{},'P3',{},'P4',{},'RCf',{},'N',{});
    
    % Plano x = 0
    AW(1).RCf = P3;
    AW(1).P1  = [0; 0; 0];
    AW(1).P2  = [0; y; 0];
    AW(1).P3  = [0; y; z];
    AW(1).P4  = [0; 0; z];
    n = cross(AW(1).P1-AW(1).P2, AW(1).P3-AW(1).P2);
    AW(1).N=n/norm(n);
    % Plano x = x
    AW(2).RCf = P3;
    AW(2).P1  = [x; 0; 0];
    AW(2).P2  = [x; y; 0];
    AW(2).P3  = [x; y; z];
    AW(2).P4  = [x; 0; z];
    n = cross(AW(2).P1-AW(2).P2, AW(2).P3-AW(2).P2);
    AW(2).N=n/norm(n);
    % Plano y = 0
    AW(3).RCf = P3;
    AW(3).P1  = [0; 0; 0];
    AW(3).P2  = [0; 0; z];
    AW(3).P3  = [x; 0; z];
    AW(3).P4  = [x; 0; 0];
    n = cross(AW(3).P1-AW(3).P2, AW(3).P3-AW(3).P2);
    AW(3).N=n/norm(n);
    % Plano y = y
    AW(4).RCf = P3;
    AW(4).P1  = [0; y; 0];
    AW(4).P2  = [0; y; z];
    AW(4).P3  = [x; y; z];
    AW(4).P4  = [x; y; 0];
    n = cross(AW(4).P1-AW(4).P2, AW(4).P3-AW(4).P2);
    AW(4).N=n/norm(n);
    % Plano z = 0
    AW(5).RCf = P3;
    AW(5).P1  = [0; 0; 0];
    AW(5).P2  = [x; 0; 0];
    AW(5).P3  = [x; y; 0];
    AW(5).P4  = [0; y; 0];
    n = cross(AW(5).P1-AW(5).P2, AW(5).P3-AW(5).P2);
    AW(5).N=n/norm(n);
    % Plano z = z
    AW(6).RCf = P3;
    AW(6).P1  = [0; 0; z];
    AW(6).P2  = [x; 0; z];
    AW(6).P3  = [x; y; z];
    AW(6).P4  = [0; y; z];
    n = cross(AW(6).P1-AW(6).P2, AW(6).P3-AW(6).P2);
    AW(6).N=n/norm(n);
else
    
    if(isempty(AW))
        AW = struct('P1',{},'P2',{},'P3',{},'P4',{},'RCf',{},'N',{});
        N = 1;
    else
        N = length(AW)+1;
    end
    
    AW(N).RCf = RCf;
    
    [M,A]=size(P1);
    if(A==1 && M==3)
        AW(N).P1=P1;
    elseif(A==1 && M==2)
        AW(N).P1=[P1;0];
    elseif(A==3 && M==1)
        AW(N).P1=P1';
    elseif(A==2 && M==1)
        AW(N).P1=[P1 0]';
    else
        error('P1 point not valid.');
    end
    
    [M,A]=size(P2);
    if(A==1 && M==3)
        AW(N).P2=P2;
    elseif(A==1 && M==2)
        AW(N).P2=[P2;0];
    elseif(A==3 && M==1)
        AW(N).P2=P2';
    elseif(A==2 && M==1)
        AW(N).P2=[P2 0]';
    else
        error('P2 point not valid.');
    end
    
    [M,A]=size(P3);
    if(A==1 && M==3)
        AW(N).P3=P3;
    elseif(A==1 && M==2)
        AW(N).P3=[P3;0];
    elseif(A==3 && M==1)
        AW(N).P3=P3';
    elseif(A==2 && M==1)
        AW(N).P3=[P3 0]';
    else
        error('P3 point not valid.');
    end
    
    [M,A]=size(P4);
    if(A==1 && M==3)
        AW(N).P4=P4;
    elseif(A==1 && M==2)
        AW(N).P4=[P4;0];
    elseif(A==3 && M==1)
        AW(N).P4=P4';
    elseif(A==2 && M==1)
        AW(N).P4=[P4 0]';
    else
        error('P4 point not valid.');
    end
    
    if(test_quad(P1,P2,P3,P4))
    elseif(test_quad(P1,P3,P2,P4))
        aux = P2;
        P2 = P3;
        P3 = aux;
    elseif(test_quad(P1,P2,P4,P3))
        aux = P3;
        P3 = P4;
        P4 = aux;
    else
        error('P1, P2, P3 and P4 not are vertices of a quadrilateral')
    end
    
    n = cross(AW(N).P1-AW(N).P2, AW(N).P3-AW(N).P2);
    
    AW(N).N=n/norm(n);
end

end

function r = test_quad(x1,x2,x3,x4)
%
% Verify if x1, x2, x3 and x4 define a quadrilateral.


t=zeros(1,4);

t(1) = acosd(dot(x1-x2,x3-x2)/(norm(x1-x2)*norm(x3-x2)));
t(2) = acosd(dot(x2-x3,x4-x3)/(norm(x2-x3)*norm(x4-x3)));
t(3) = acosd(dot(x3-x4,x1-x4)/(norm(x3-x4)*norm(x1-x4)));
t(4) = acosd(dot(x4-x1,x2-x1)/(norm(x4-x1)*norm(x2-x1)));

r = round(sum(t)) == 360;
end
