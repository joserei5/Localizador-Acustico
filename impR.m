function I = impR(file, R, C, varargin)
%
%  impR - Determines the impulse and frequency response of each real source.
%
%  I = impR(file, R, C)
%
%  INPUT:
%    file  - the name of the file with the room data;
%    R     - receiver data (created with the function addReceiver);
%    C     - Channel object
%
%
%  Copyright 2012 LocUS Project.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 2.2$    $Date: 2021/01/20$

if(nargin~=3)
    error('Not enough input arguments.');
end

load(file,'R_ind','AW','S_coord','S_orig','S_Ex','S_Ey','S_Ez','A_TF','fs','mR','flim','A');


if(~isempty(flim))
    C.fl = flim(1);
    C.fh = flim(2);
end

C.fs = fs;

% Sources validation
M_Sv = valid_S(R_ind,AW,S_coord,S_orig,R.C,mR);
iSv = logical(M_Sv);
pSv = find(iSv);

if(sum(M_Sv)>0)
    
    RSc = S_coord(:,iSv);
    
    % Valid sources index
    Sv = find(M_Sv);
    
    % Original source TF index
    iTF = ceil(Sv.*sum(R_ind==0)./length(R_ind));
    
    
    % Distance
    d = calc_d(RSc,R.C);
        
    v_SR = zeros(size(RSc));
    
    v_SR(1,:) = R.C(1) -RSc(1,:);
    v_SR(2,:) = R.C(2) -RSc(2,:);
    v_SR(3,:) = R.C(3) -RSc(3,:);
    
    RS_Ex = S_Ex(:,iSv);
    RS_Ey = S_Ey(:,iSv);
    RS_Ez = S_Ez(:,iSv);
    
    n_V_SR_x = dot(v_SR,RS_Ex,1);
    n_V_SR_y = dot(v_SR,RS_Ey,1);
    n_V_SR_z = dot(v_SR,RS_Ez,1);
    
    ATT = struct('h',{},'p',{});
    
    v_RS = -v_SR;
    
    mh = 0; % max value of h
    mp = 0; % max value of p
    
    for n=1:length(n_V_SR_x)
        
        % Wall filter
        ATT(n).h = A(pSv(n)).h;
        ATT(n).p = A(pSv(n)).p;
        
        % Source filter
        [theta,phi,r]=cart2sph(n_V_SR_x(n),n_V_SR_y(n),n_V_SR_z(n));
        
        if(isempty(flim))
            [h, p] = A_TF{iTF(n)}(theta,phi,fs);
        else
            [h, p] = A_TF{iTF(n)}(theta,phi,fs,flim(1),flim(2));
        end
        
        ATT(n).h = conv(ATT(n).h,h);
        ATT(n).p = ATT(n).p+p;
        
        % Receiver filter
        [theta,phi,r]=cart2sph(dot(v_RS(:,n),R.D(:,1)),dot(v_RS(:,n),R.D(:,2)),dot(v_RS(:,n),R.D(:,3)));
        
        if(isempty(flim))
            [h, p] = R.Tf(theta,phi,fs);
        else
            [h, p] = R.Tf(theta,phi,fs,flim(1),flim(2));
        end
        
        ATT(n).h = conv(ATT(n).h,h);
        ATT(n).p = ATT(n).p+p;
        
        % Air Filter
        [h, p] = impulse(C,d(n));
        
        ATT(n).h = conv(ATT(n).h,h);
        ATT(n).p = ATT(n).p+p;
        
        if(length(ATT(n).h)>mh)
            mh = length(ATT(n).h);
        end
        
        if(ATT(n).p > mp)
            mp = ATT(n).p;
        end
        
    end
    
%     I = zeros(max(iTF),mp+mh);
    I = zeros(sum(R_ind==0),mp+mh);


    % Impulse response for each source
    for n=1:length(n_V_SR_x)
        Nh = length(ATT(n).h);
        I(iTF(n),ATT(n).p+1:ATT(n).p+Nh) = I(iTF(n),ATT(n).p+1:ATT(n).p+Nh)+ATT(n).h;
    end
    
    
else
    I=[];
end
end

function M_Sv = valid_S(R_ind,Awall,S_coord,S_orig,R_coord,MAX_RI)
%
% valid_S validate all the virtual sources. If the path of some ray cross
% some object or surface this ray is not valid so the virtual source do not
% exist.

Nsource = size(R_ind,2);
Nwall = size(Awall,2);
M_Sv = ones(1,Nsource);

% h = waitbar(0,'Virtual sources validation...');

NBAR = (MAX_RI*Nsource+Nsource);
jump_bar = round(NBAR/10);
nbar = 0;

for T_RI = 1:MAX_RI
    ind = Nsource;
    while ind>0
        nbar = nbar + 1;
        if(rem(nbar,jump_bar)==0)
%             waitbar(nbar/NBAR,h);
        end
        if(R_ind(ind)==T_RI)
            ind2 = ind;
            ind22 = ind;
            T_Pr = R_coord;
            while(R_ind(ind22)>0 && M_Sv(ind)) %
                if(R_ind(ind2)==R_ind(ind22)-1)
                    if(S_orig(ind22) ~= 0)
                        T_Pr_old = T_Pr;
                        [T_Pr M_Sv_aux1] = inter_PP(T_Pr,S_coord(:,ind22),Awall(S_orig(ind22)).N,...
                            Awall(S_orig(ind22)).P1); % Test if the ray interset the plan that originate the reflection
                        if(M_Sv_aux1)
                            for ind_s=1:Nwall; % test if the ray intersect another surface.
                                if(ind_s~= S_orig(ind22))
                                    [T_Prx, M_Sv_aux2]= inter_PP(T_Pr_old,T_Pr,Awall(ind_s).N,Awall(ind_s).P1);
                                    if(M_Sv_aux2)
                                        M_Sv(ind) = M_Sv(ind) & ~inquad(T_Prx,Awall(ind_s).P1,Awall(ind_s).P2,...
                                            Awall(ind_s).P3,Awall(ind_s).P4);
                                    end
                                end
                            end
                            
                            M_Sv(ind) = M_Sv(ind) & inquad(T_Pr,Awall(S_orig(ind22)).P1,...
                                Awall(S_orig(ind22)).P2,Awall(S_orig(ind22)).P3,Awall(S_orig(ind22)).P4);
                            
                            % %%%%%%%%%%%%%%
                            if(R_ind(ind22)==1 && M_Sv(ind)) % for indice 1 test also if the ray to the source intersect any surface
                                
                                for ind_s=1:Nwall % test if the ray intersect another surface.
                                    if(ind_s~= S_orig(ind22))
                                        [T_Prx, M_Sv_aux2]= inter_PP(T_Pr,S_coord(:,ind2),Awall(ind_s).N,Awall(ind_s).P1);
                                        if(M_Sv_aux2 || sum(abs(T_Pr-T_Prx))<eps)
                                            %keyboard
                                            M_Sv(ind) = M_Sv(ind) & ~inquad(T_Prx,Awall(ind_s).P1,Awall(ind_s).P2,...
                                                Awall(ind_s).P3,Awall(ind_s).P4);
                                        end
                                    end
                                end
                            end
                            % %%%%%%%%%%%%%%
                            
                        else
                            M_Sv(ind)=false;
                        end
                    else
                        keyboard
                    end
                    ind22=ind2;
                end
                ind2=ind2-1;
            end
        end
        ind = ind -1;
        
    end
end

for ind=1:Nsource
    nbar = nbar + 1;
    if(rem(nbar,jump_bar)==0)
%        waitbar(nbar/NBAR,h);
    end
    if(R_ind(ind)==0)
        for ind_s=1:Nwall; % test if the ray intersect another surface.
            T_Prx = inter_PP(R_coord,S_coord(:,ind),Awall(ind_s).N,Awall(ind_s).P1);
            if(inseg(R_coord,S_coord(:,ind),T_Prx))
                M_Sv(ind) = M_Sv(ind) & ~inquad(T_Prx,Awall(ind_s).P1,...
                    Awall(ind_s).P2,Awall(ind_s).P3,Awall(ind_s).P4);
            end
        end
    end
end
% close(h);
end

function [Pi, Pv] = inter_PP(R1,R2,N,P)
%
%  Test if the line segment intersects the surface defines by the normal
%  vector N and the point P

% a vector in the line
v = R2-R1;
%a = dot(v,N);
a = v'*N;

% equation of the surface -> X.nx+Y.ny+Z.nz=d
%d = dot(N,P);
d = N'*P;

if(a == 0) % not apply the equation
    if(d == 0) % surface in the origin
        Pi = R1;
        Pv = true;
    else
        Pi = zeros(size(R1));
        Pv = false; % not intersect
    end
else
    %k = (d-dot(R1,N))/a; % distance to the intersection.
    k = (d-R1'*N)/a;
    Pi = v*k + R1;
    Pv = inseg(R1,R2,Pi); % verify if the point is inside the line segment
end

end

function IN = inquad(P,P1,P2,P3,P4)
%
% Verify if the point P is inside the quadrilateral defines by the points
% P1, P2, P3, and P4.

A(1) = t_area(P1,P2,P);
A(2) = t_area(P2,P3,P);
A(3) = t_area(P3,P4,P);
A(4) = t_area(P4,P1,P);

AT(1) = t_area(P1,P2,P3);
AT(2) = t_area(P1,P3,P4);

IN = abs((sum(AT)-sum(A))/sum(A))<1e-3;

end

function IN = inseg(P1,P2,P)
%
%  Verify if the point P is in the line segment

IN = true;
for c=1:3
    IN = IN & ((P(c)>= P1(c) & P(c)<= P2(c) | P(c)<= P1(c) & P(c)>= P2(c)));
end
IN = IN & (P(1)~=P1(1) | P(2)~=P1(2) | P(3)~=P1(3)); % limits are now
%valid
IN = IN & (P(1)~=P2(1) | P(2)~=P2(2) | P(3)~=P2(3));
end

function A = t_area(x1,x2,x3)
%
%  Calculate the triangle area.

a = norm(x1-x2);
b = norm(x1-x3);
c = norm(x2-x3);

p = (a + b + c) / 2;

A = sqrt(p * (p - a) * (p - b) * (p - c));

end

function d=calc_d(S,R)
%
%  Determines the distance of each source to the receiver.

N=size(S,2);
d = zeros(1,N);

for(ind=1:N)
    d(ind) = norm(S(:,ind)-R);
end
end