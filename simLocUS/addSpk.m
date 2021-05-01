function S=addSpk(C,Tf,D)
%
%  addSpeaker - Add a "Speaker" to the Speaker array (S).
%
%  R = addSpk(C)
%
%   or
%
%  S = addSpk(C,Tf,D)
%
%   INPUT:
%   C     - Speaker coordinates (vector with the coordinates [x y] or [x y z]);
%   Tf    - (optional, the default is omnidirectional without transfer function)
%           transfer function handle of the transducer (impulse response)
%           if you set this input you must set the Radiation direction (D).
%   D     - Radiation direction of the source [Dx Dy] or [Dx Dy Dz] (optional,
%           the default is omnidirectional without transfer function);
%
%
%  Copyright 2020.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 1.1$    $Date: 2020/10/30$

if(nargin<1)
    error('Not enough input arguments.');
end

E = numel(C);
if(E==2)
    S.C = [C(:); 0];
elseif(E==3)
    S.C = C(:);
else
    error('Invalid receiver coordinates.');
end

if(nargin==3)
    if isa(Tf, 'function_handle')
        S.Tf = Tf;
    else
        error('Transfer function not valid.');
    end
    
    [M,A]=size(D);
    E = numel(D);
    if(E==2)
        D = [D(:); 0];
        D = orthoNVec(D);
    elseif(E==3)
        D = D(:);
        D = orthoNVec(D);
    elseif(A==3 && M==2)
        D = [D;zeros(1,3)];
    elseif(A==2 && M==3)
        D = [D';zeros(1,3)];
    elseif(E==9)
    else
        error('Invalid radiation direction of the source.');
    end
    
    D(:,1)=D(:,1)./norm(D(:,1));
    D(:,2)=D(:,2)./norm(D(:,2));
    D(:,3)=D(:,3)./norm(D(:,3));

    S.D = D;
elseif(nargin==1)
    S.Tf = @omni;
    S.D = orthoNVec([1;0;0]);
else
    error('Not enough input arguments.');
end

end

function [h,p] = omni(th,ph,fs,fl,fh)
    h = 1;
    p = 0;
end