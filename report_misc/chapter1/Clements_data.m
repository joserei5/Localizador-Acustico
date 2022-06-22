clear;clc;close all;

%% Data1
ae.aegypti.f_T = [    linspace(508-66,508+66,5)'; ...
                    linspace(414,557,5)';       ...
                    linspace(458,470,5)';       ];

ae.aegypti.T = [22 23 25];
g1 = repmat({'22'},5,1);
g2 = repmat({'23'},5,1);
g3 = repmat({'25'},5,1);
g = [g1; g2; g3];
           
%% Figures1
f1=figure();
% sbplt1 = subplot(1,2,1,'Parent',f1);
boxplot(ae.aegypti.f_T, g)
ylabel('Frequency (Hz)')
xlabel('Temperature (Celsius)')
title('Aedes aegypti wing-beat frequency')

%% Data2
ae.aegypti.f = round(linspace(414,508+66,10))';
ae.albopictus.f = round(linspace(536,544,10))';
ae.diantaeus.f = round(linspace(317,434,10))';
ae.communis.f = round(linspace(350-9,350+9,10))';
ae.punctor.f = round(linspace(308-10,308+10,10))';
ae.triseriatus.f = round(linspace(388-33,388+33,10))';
ae.vexans.f = round(linspace(300,350,10))';

ae_all_sp = [   ae.aegypti.f;       ...
                ae.albopictus.f;    ...
                ae.diantaeus.f;     ...
                ae.communis.f;      ...
                ae.punctor.f;       ...
                ae.triseriatus.f;     ...
                ae.vexans.f;          ];

g1 = repmat({'aegypti'},10,1);
g2 = repmat({'albopictus'},10,1);
g3 = repmat({'diantaeus'},10,1);
g4 = repmat({'communis'},10,1);
g5 = repmat({'punctor'},10,1);
g6 = repmat({'triseriatus'},10,1);
g7 = repmat({'vexans'},10,1);
g = [g1; g2; g3; g4; g5; g6; g7];

%% Figures2
% figure(f1);
f2 = figure();
% sbplt2 = subplot(1,2,2,'Parent',f1);
boxplot(ae_all_sp, g, 'orientation', 'horizontal')
ylabel('Species')
xlabel('Frequency (Hz)')
title('Aedes wing-beat frequency')

% % adjust subplot widths
% width_adjust = 0.05;
% sbplt1.Position(3) = sbplt1.Position(3) - width_adjust;
% sbplt2.Position(1) = sbplt2.Position(1) - width_adjust;
% sbplt2.Position(3) = sbplt2.Position(3) + 0.1;

%% Data3
an.earlei.f = round(linspace(190,250,10))';
an.subpictus.f = round(linspace(330,385,10))';
an.arabiensis.f = round(linspace(360,520,10))';
an.gambiae.f = round(linspace(420,600,10))';
an.melas.f = round(linspace(368,377,10))';
an.merus.f = round(linspace(372,288,10))';

an_all_sp = [   an.earlei.f;        ...
                an.subpictus.f;     ...
                an.arabiensis.f;	...
                an.gambiae.f;       ...
                an.melas.f;         ...
                an.merus.f;         ];

g1 = repmat({'earlei'},10,1);
g2 = repmat({'subpictus'},10,1);
g3 = repmat({'arabiensis'},10,1);
g4 = repmat({'gambiae'},10,1);
g5 = repmat({'melas'},10,1);
g6 = repmat({'merus'},10,1);
g = [g1; g2; g3; g4; g5; g6];

%% Figures3
f3=figure();
% sbplt1 = subplot(1,2,1,'Parent',f1);
boxplot(an_all_sp, g, 'orientation', 'horizontal')
ylabel('Species')
xlabel('Frequency (Hz)')
title('Anopheles wing-beat frequency')


%% Change export renderer
set(f1,'Renderer','painters')
set(f2,'Renderer','painters')
set(f3,'Renderer','painters')