function AM = addMic(AM,C,Tf,D)
%
%  addMic - Add a "Microphone" to the Microphone array (AM).
%
%  AM = addMic(AS, C)
%
%   or
%
%  AM = addMic(AS, C, Tf, D)
%
%   INPUT:
%   AM    - Mic. array that you want to add the microphone;
%   C     - Mic. coordinates (vector with the coordinates [x y] or [x y z]);
%   Tf    - (optional, the default is omnidirectional without transfer function)
%           ransfer function handle of the transducer (impulse response)
%           if you set this input you must set the Radiation direction (D).
%   D     - Radiation direction of the source [Dx Dy] or [Dx Dy Dz] (optional,
%           the default is omnidirectional without transfer function);
%
%
%  Copyright 2020.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 1.1$    $Date: 2020/10/30$

if(nargin<2)
    error('Not enough input arguments.');
end

if(isempty(AM))
    AM = struct('C',{},'Tf',{},'D',{},'h',{},'p',{});
    % C - cordinates
    % Tf - Transferfunction
    % D - direction
    % h - filter coeficients due to reflections (1 -> direct path)
    % p - filter position (0 -> direct path)
    
    N = 1;
else
    N = length(AM)+1;
end

if(nargin==4)
    if isa(Tf, 'function_handle')
        AM(N).Tf=Tf;
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
    
    AM(N).D=D;
    
elseif(nargin==2)
    AM(N).Tf = @omni;
    AM(N).D = orthoNVec([1;0;0]);
else
    error('Not enough input arguments.');
end


[M,A]=size(C);
if(A==1 && M==3)
    AM(N).C = C;
elseif(A==1 && M==2)
    AM(N).C = [C;0];
elseif(A==3 && M==1)
    AM(N).C = C';
elseif(A==2 && M==1)
    AM(N).C = [C 0]';
else
    error('Coordinate vector not valid');
end

AM(N).h=1;
AM(N).p=0;
end

function [h,p] = omni(th,ph,fs,fl,fh)
    h = 1;
    p = 0;
end
