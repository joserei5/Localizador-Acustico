function [f,a]=displayRoom(Ffile,op)
%
%  displayRoom - Show the surfaces and sources in 3D.
%
%  [f,a] = displayRoom(Ffile);
%  [f,a] = displayRoom(Ffile,'HideVS');
%
%  INPUT:
%   Ffile - File created with the function makeFile;
%   HideVS - Hide the Virtual Sources
%
%  OUTPUT:
%   f - Figure handle;
%   a - Axis handle;
%
%  Copyright 2012 LocUS Project.
%  Written by:  Daniel Filipe Albuquerque.
%  $Revision: 2.3$    $Date: 2012/05/09$

    load (Ffile,'R_ind','S_coord','AW','S_Ex','mR');

    f=figure;
    a=axes;
    
    colormap(hot)
    C = colormap;
    NColor = size(C,1)-15;

    if(~exist('op','var'))
        op = 'none';
    end
    
    if(strcmp(op,'HideVS'))
        mR=0;
    end
    
     for r=0:mR;
        plot3(S_coord(1,R_ind==r),S_coord(2,R_ind==r),S_coord(3,R_ind==r)...
        ,'o','markersize',18-12*r/(mR+1),'markeredgecolor',[0 0 0],'markerfacecolor',...
        C(round(NColor*(r+1)/(mR+1)),:))
        hold on
     end
     grid on;
     
     
     
     for r=1:size(S_coord,2);
         if(strcmp(op,'HideVS'))
             if(R_ind(r)==0)
                 line([S_coord(1,r) S_coord(1,r)+S_Ex(1,r)],[S_coord(2,r) S_coord(2,r)+S_Ex(2,r)],...
                     [S_coord(3,r) S_coord(3,r)+S_Ex(3,r)],'Color',C(round(NColor*(R_ind(r)+1)/(mR+1)),:),...
                     'LineWidth',2,'Marker','.','MarkerSize',16,'MarkerFaceColor',C(round(NColor*(R_ind(r)+1)/(mR+1)),:));
             end
             
         else
             
             
             line([S_coord(1,r) S_coord(1,r)+S_Ex(1,r)],[S_coord(2,r) S_coord(2,r)+S_Ex(2,r)],...
                 [S_coord(3,r) S_coord(3,r)+S_Ex(3,r)],'Color',C(round(NColor*(R_ind(r)+1)/(mR+1)),:),...
                 'LineWidth',2,'Marker','.','MarkerSize',16,'MarkerFaceColor',C(round(NColor*(R_ind(r)+1)/(mR+1)),:));
         end
     end
     
     N_Sur = size(AW,2);
     
     colormap(hsv)
    C = colormap;
    NColor = size(C,1);
    
    for n = 1:N_Sur;
        fill3([AW(n).P1(1) AW(n).P2(1) AW(n).P3(1) AW(n).P4(1)],...
              [AW(n).P1(2) AW(n).P2(2) AW(n).P3(2) AW(n).P4(2)],...
              [AW(n).P1(3) AW(n).P2(3) AW(n).P3(3) AW(n).P4(3)],...
              C(round(NColor*n/N_Sur),:),'FaceAlpha',0.4,...
              'EdgeColor',[1,1,1],'EdgeAlpha',0.6)
    end
        
        xlabel('x','FontSize',14,'FontName','Courier New');
        ylabel('y','FontSize',14,'FontName','Courier New');
        zlabel('z','FontSize',14,'FontName','Courier New');
        set(a,'FontSize',14,'FontName','Courier New','Box','off','FontAngle','italic');
        axis equal
end