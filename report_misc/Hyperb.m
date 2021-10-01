clear;clc;close all;

%%
syms x y z

% Microphone (focus) spacing in centimeters
Delta_x = 29.5 * 1e-2; 
% Microphone (focus) x-axis location
c = Delta_x * 1/2;
% Hyperbole Vertices distance
const = 0:c/4:Delta_x;
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
figure;
    fimplicit(hyperb)
    xlim([-c-.25*c c+.25*c]);
    ylim([-1 1])
    title(  '$\sqrt{(x+c)^2+y^2} - \sqrt{(x-c)^2+y^2} = \pm \: 2a$',...
            '$\Leftrightarrow (x/a)^2 - (y/b)^2 == 1 \quad \wedge \quad b^2=c^2-a^2$',...
            'Interpreter','latex')
    legend(cellstr(num2str(const', 'a=%.4f')))
    
figure
    fimplicit([exp==const exp==-const])
    xlim([-c-.25*c c+.25*c]);
    ylim([-1 1])
    title(  '$\sqrt{(x+c)^2+y^2} - \sqrt{(x-c)^2+y^2} = \pm \: 2a$',...
            '$\Leftrightarrow (x/a)^2 - (y/b)^2 == 1 \quad \wedge \quad b^2=c^2-a^2$',...
            'Interpreter','latex')