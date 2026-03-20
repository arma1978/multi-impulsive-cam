clearvars; close all; clc; format longG;
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultUicontrolFontName','Times', 'DefaultUicontrolFontSize', 14);
set(0,'DefaultUitableFontName','Times', 'DefaultUitableFontSize', 14);
set(0,'DefaultTextFontName','Times', 'DefaultTextFontSize', 14);
set(0,'DefaultUipanelFontName','Times', 'DefaultUipanelFontSize', 14);

projectRoot = fileparts(fileparts(mfilename('fullpath')));
cd(projectRoot);
addpath(fullfile(projectRoot,'src','routines'));
addpath(fullfile(projectRoot,'src','io'));
addpath(fullfile(projectRoot,'src','frames'));
addpath(fullfile(projectRoot,'utils'));
addpath(fullfile(getenv('HOME'),'Dropbox/Work/mosek/11.1/toolbox/r2022b'));

name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_Kepmiss_distance.mat');
load(name);
outKep = out; 
name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_J2miss_distance.mat');
load(name);
outJ2 = out; 

A = [outJ2.ind];
B = [outKep.ind];
[c,ia,ib] = intersect(A,B); 
outJ2 = outJ2(ia);
outKep = outKep(ib);
%% plotting
%% plotting the results
res = '-r500';
type = '-dpng';
path = fullfile(projectRoot,'output','Figures');
if ~exist(path,'dir'), mkdir(path); end

lowerPer = 95;
upperPer = 100;
disc = 100;

diffdvKepJ2 = abs([outKep.dv]-[outJ2.dv]);

[a,b] = max(diffdvKepJ2)

figure(); hold on; grid on; box on
nameFig = 'diffdvKepJ2';
histogram(rmoutliers(diffdvKepJ2*km^2,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$||\Delta {\bf v}_{J_2-J_4} - \Delta {\bf v}_{\textrm{kep}}||$ [mm/s]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

bplaneKep = [outKep.rrb2dRef];
bplaneJ2 = [outJ2.rrb2dRef];
diffbplaneKepJ2 = bplaneKep-bplaneJ2; 
diffbplaneKepJ2 = sqrt(sum(diffbplaneKepJ2.*diffbplaneKepJ2,1))*km;

figure(); hold on; grid on; box on
nameFig = 'diffbplaneKepJ2';
histogram(rmoutliers(diffbplaneKepJ2,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$||\Delta {\bf r}_{CA,J_2-J_4} - \Delta {\bf r}_{CA,\textrm{kep}}||$ [m]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

bplaneLin = [outJ2.rrb2dLin];
bplaneJ2 = [outJ2.rrb2dRef];
diffbplaneLin = bplaneLin-bplaneJ2; 
diffbplaneLin = sqrt(sum(diffbplaneLin.*diffbplaneLin,1))*km;

figure(); hold on; grid on; box on
nameFig = 'diffbplaneLin';
histogram(rmoutliers(diffbplaneLin,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$||\Delta {\bf r}_{CA} - \Delta {\bf r}_{CA,\textrm{lin}}||$ [m]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

diffdvLin = abs([outJ2.dvLin]-[outJ2.dv]);

figure(); hold on; grid on; box on
nameFig = 'diffdvLin';
histogram(rmoutliers(diffdvLin*km^2,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$||\Delta {\bf v} - \Delta {\bf v}_{\textrm{lin}}||$ [mm/s]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);


