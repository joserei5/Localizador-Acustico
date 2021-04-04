function makeFile(Name, AW, AS, mR, fs, fl, fh)
%
%  makeFile - Create a file with all the important Data of the system for 
%             in the simulator.
%
%  makeFile(Name, AW, AS, ind_r)
%
%  INPUT:
%   Name    - Name of the file;

%   AW      - Array with all the surfaces;
%   AS      - Array with all the sources;
%   mR      - Maximum number of reflections to consider.
%   fs      - Sample frequency to be used (It is recomended to use 4 times 
%             more than the maximum signal frequency).
%   fl-fh   - (optional) frequency bounds for frequency response optimization
%             (default) 0-fs/2;
%
%  Copyright 2011 LocUS Project.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 2.0$    $Date: 2011/11/18 17:50$

    if(nargin==5)
        flim = [];
    elseif(nargin==7)
        if(fl>fh)
            error('Invalid frequency bounds.')
        end
        flim    = [fl, fh];
    else
        error('Not enough input arguments.');
    end

    if (~isempty(AW))
        ATREE = calc_T(AS, AW, mR, fs, flim);
    else
        ATREE = calc_T(AS, AW,  0, fs, flim);
    end
    
    A = c_ASource(ATREE);
    
    A_TF = c_A_TF(ATREE);

    %%%%%%%%%%%%%%%%% conversion %%%%%%%%%%%%%%%%%%%%%
    R_ind   = [A.R];
    S_orig  = [A.S];
    S_coord = [A.C];
    aux     = [A.D];
    S_Ex    = aux(:,1:3:end);
    S_Ey    = aux(:,2:3:end);
    S_Ez    = aux(:,3:3:end);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(Name,'R_ind','S_orig','S_coord','S_Ex','S_Ey','S_Ez','ATREE','AW','AS','A_TF','A','mR','fs','flim');

end

function C = c_ASource(A)
%
%  Read the tree and put the data into a matrix.


    C = struct('R',{},'S',{},'C',{},'D',{},'h',{},'p',{});
    % R - reflection index
    % S - wall that creat the reflection (0 -> direct path)
    % C - coordinates
    % D - direction
    % h - filter coeficient
    % p - filter position

    N = size(A,2);
 
    for n=1:N
        p = length(C)+1;
        C(p).R = 0;
        C(p).S = A(n).S_orig;
        C(p).C = A(n).Source.C;
        C(p).D = A(n).Source.D;
        C(p).h = A(n).Source.h;
        C(p).p = A(n).Source.p;
        C = coord_I(C,A(n),0);
    end
end

function C=coord_I(C,A,R)
%
%  Recursive function for use in the function c_ASource.


    N = size(A.Imag,2);
    if(N>0)
        for n=1:N
            p = length(C)+1;
            C(p).R = R+1;
            C(p).S = A.Imag(n).S_orig;
            C(p).C = A.Imag(n).Source.C;
            C(p).D = A.Imag(n).Source.D;
            C(p).h = A.Imag(n).Source.h;
            C(p).p = A.Imag(n).Source.p;
            C=coord_I(C,A.Imag(n),R+1);
        end
    end
end

function A_TF = c_A_TF(ATREE)
%
% Create an array with the transfer functions of all real sources.

    N = size(ATREE,2);
    A_TF = cell(1,N);
    for n=1:N
        A_TF{n} =ATREE(n).Source.Tf;
    end
end

function A = calc_T(ASource,AWall,R,fs,flim)
%
%  Compute the tree up to R reflection coefficiente.

    Mf = size(ASource,2);
    A=[];
    for mf=1:Mf;
        A=[A calcI(ASource(mf),0,AWall,R,fs,flim)];
    end
end

function A=calcI(F,S,AWall,R,fs,flim)
%
%  Recursive function for use in the function calc_T.

    Ms = size(AWall,2);
    A = struct('Source',[],'S_orig',[],'Imag',[]);
    A.Source = F;
    A.S_orig = S;
    
    if(R>0)    
        N = [AWall.N];
        P = [AWall.P1];
        E = calc_imag(F.C,N,P);
        
        I = struct('C',{},'Tf',{},'D',{},'h',{},'p',{});
        for ms=1:Ms
            if(ms~=S)
                I(ms) = F;
                I(ms).C = E(:,ms);
                I(ms).D = F.D+2*N(:,ms)*(-N(:,ms)'*F.D);
                
                if(isa(AWall(ms).RCf,'function_handle'))
                    if(isempty(flim))
                        [h,p] = AWall(ms).RCf(fs);
                        
                    else
                        [h,p] = AWall(ms).RCf(fs,flim(1),flim(2));
                    end
                else
                    h = AWall(ms).RCf;
                    p = 0;
                end
                I(ms).h = conv(F.h,h);
                I(ms).p = F.p+p;
                
                A.Imag = [A.Imag calcI(I(ms),ms,AWall,R-1,fs,flim)];
            end
        end
    end
    
    
end

function I=calc_imag(S,N,P)
%
%  Find the virtual image of source S on the surface defines for the
%  normal vector N and the point P.

    Mw = size(N,2);
    Ms = size(S,2);
    
    I = zeros(3,Mw,Ms);
    
    for ms=1:Ms;
        AUX_S = S(:,ms)*ones(1,Mw);
        I(:,:,ms)=AUX_S+2*N*diag(diag(N'*(P-AUX_S)));
    end
end
