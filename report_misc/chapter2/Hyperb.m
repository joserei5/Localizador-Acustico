clear;clc;close all;

%% GLOBAL KEYS
% activate alternative expression
Talt = 0;
% activate branch demo
Tbd = 0;

%%
syms x y z

% Microphone (focus) spacing in centimeters
Delta_x = 29.5 * 1e-2; 
% Microphone (focus) x-axis location
c = Delta_x * 1/2;
% Hyperbole Vertices distance
const = 0:Delta_x/10:Delta_x;
% Hyperbole Vertex distance to origin
a = const/2;
% Hyperbole distance expression
exp = sqrt((x+c)^2 + y^2)-sqrt((x-c)^2 + y^2);
% Hyperbole alternative expression
b = sqrt(c^2 - a.^2);
hyperb = (x./a).^2 - (y./b).^2 == 1;


% VER 2 EXPRESSOES
% NAO APARECEM AS LINHAS HORIZ E VERT NA ULTIMA

%%
if Talt
    figure;
        ax = fimplicit(hyperb);
        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        xlim([-c-.25*c c+.25*c]);
        ylim([-1 1])
    %     title(  '$\sqrt{(x+c)^2+y^2} - \sqrt{(x-c)^2+y^2} = \pm \: 2a$',...
    %             '$\Leftrightarrow (x/a)^2 - (y/b)^2 == 1 \quad \wedge \quad b^2=c^2-a^2$',...
    %             'Interpreter','latex')
        legend(cellstr(num2str(transpose(b)*1e2, '%.2f')))
end
    
figure
    ax = fimplicit(abs(exp)==const);
    xlabel('x')
    ylabel('y')
    xlim([-c-.25*c c+.25*c]);
    ylim([-1 1])
%     title(  '$\sqrt{(x+c)^2+y^2} - \sqrt{(x-c)^2+y^2} = \pm \: 2a$',...
%             '$\Leftrightarrow (x/a)^2 - (y/b)^2 == 1 \quad \wedge \quad b^2=c^2-a^2$',...
%             'Interpreter','latex')
    legend(cellstr(num2str(transpose(const)*1e2, '%.2f')))
    
%%
if Tbd
    curveN = 3;
    figure;
        ax = fimplicit([exp==const(curveN) exp==-const(curveN)],'k');
        hold on;
        plot(c,0,'b+');
        plot(-c,0,'b+');
        hold off;
        xlabel('x')
        ylabel('y')
        xlim([-1 1]);
        ylim([-1 1])
end