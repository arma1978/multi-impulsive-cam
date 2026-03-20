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
runtimeDir = fullfile(projectRoot,'runtime');
outputDir = fullfile(projectRoot,'output');
if ~exist(runtimeDir,'dir'), mkdir(runtimeDir); end

%% DATA %% all units are km, s,
load(fullfile(projectRoot,'input','data.mat'))
km = 1000;
param.mu = 398600.4418;
param.rE = 6378.137;
param.J2 = 1.08262668e-3;
param.Lsc = param.rE;
param.Vsc = sqrt(param.mu/param.Lsc);
param.Tsc = param.Lsc/param.Vsc;
param.disc = 301;
mass = 300; % kg
Tmax = 30e-3; %1000e-3; %15e-3; % 50e-3; %1000e-3;%13.0007e-3; %millinew

%% user input data
param.prop = 'Kep';  % J2
param.nOrbits = 2; % time to maneuver in terms of period of sc
param.ind = 889;     % which test case
param.order = 2;
param.nvar  = 7;
param.dt = 60; % in seconds
param.maxImpNum = 170; % maximum number of impulses
param.dvMax = Tmax /mass*param.dt/1000;     % single dv km/s
PcLim = 1e-6;
maxPcLim = 1e-4;
minDistLim = 2;
sim = 1;

switch sim
    case 1
        tag = 'miss_distance';
        param.minDist = minDistLim; % rivedere!! provide this for min distance
        param.PcMax = maxPcLim;
        param.Pc = [];
    case 2
        tag = 'risk';
        param.Pc = PcLim;
        param.minDist = [];
        param.PcMax = [];
    case 3
        tag = 'max_risk';
        param.PcMax = maxPcLim;
        param.minDist = [];
        param.Pc = [];
end

%% parameters setting
param.outfile = fullfile(outputDir,'ResultsSCVX','sol');
param.fname   = fullfile(runtimeDir,'maps.dat');
param.projectRoot = projectRoot;
param.runtimeDir = runtimeDir;
if ~exist(fileparts(param.outfile),'dir'), mkdir(fileparts(param.outfile)); end

%% extract the test case and prepare the file for multi-impulse.
param.data = data;
% param.data = sortrows(data,'v^* [km/s]', 'ascend');
param = selectTestCaseFiltered(param);

%% run the main code (system for windows)
if strcmp(param.prop,'J2')
   % !./runPropMultiMapsFullPKep
    status = system(fullfile(projectRoot,'bin','statePropMultiMapsKep'));
    assert(status==0,'statePropMultiMapsKep execution failed.');
    dum = load(fullfile(runtimeDir,'xs0Kep.dat'));
    param.xx0s = dum(1:6);
else
   % !./runPropMultiMapsFullP
    status = system(fullfile(projectRoot,'bin','statePropMultiMapsFullP'));
    assert(status==0,'statePropMultiMapsFullP execution failed.');
    dum = load(fullfile(runtimeDir,'xs0J2.dat'));
    param.xx0s = dum(1:6);
end
%% load the output
N = min(floor(param.t2TCA/param.dt),param.maxImpNum);
nmaps = (N+2)*6+2+4+4; % N*6 maps, 6*2 maps at encounter, tca, mhsqr, maha elements, cov
DAx = loadPoly(param.fname,param.nvar,nmaps,0);

%% prepare for optimisation
nOptVarDvComp = N*3;
nOptVarDvMag  = N;
nOptVarState  = 2;
nOptVar = nOptVarDvComp + nOptVarDvMag + nOptVarState;
param.nOptVar = nOptVar;
param.nOptVarDvComp = nOptVarDvComp;
param.nOptVarDvMag  = nOptVarDvMag;
param.nOptVarState  = nOptVarState;
param.ndv = nOptVarDvComp/3;

%% pre-calculating relevant quantities
param = refTraj(DAx,param);
param = linearMaps(DAx,param);
param = ellipseTargets(param);

[sqrMaha0,tca0,dv0,dvv0,xx0,...
    xxNoMan0,xxsEnc0,xxdEnc0,rrb2d0,...
    Pb_2D0, maxPc0, Pc0] = validate(zeros(nOptVarDvComp,1),DAx,param);

%% convex solution 1st round _Lin
[outLin,xxOptLin,param] = linConvexSolveEllipse(param);

%% scatter plot! 
A =[ outLin(:).xxOnEllipse];
B = [outLin(:).dvTot]*km^2;
scatter([A(2,:)],[A(1,:)],[] ,B,'filled');
hold on 
%plot([A(2,116)],[A(1,116)], 'ws')
%plot([A(2,266)],[A(1,266)], 'ks')
[~,maxInd] = findpeaks(B);
[minDv,minInd] = findpeaks(-B);
[a,b] = max(minDv);
%plot(A(2,maxInd),A(1,maxInd), 'k^','markersize',10, 'MarkerFaceColor', 'w' );
plot(A(2,minInd),A(1,minInd), 'kp','markersize',10, 'MarkerFaceColor', 'w');
plot(A(2,minInd(b)),A(1,minInd(b)), 'kh','markersize',10, 'MarkerFaceColor', 'w');
ylabel('$\xi$ [km]')
xlabel('$\zeta$ [km]')
colorbar
title('$\Delta v$ [mm/s]')
caxis([min(B) max(B)])
axis equal
grid on 
box on 


% %% save workspace and output file
% clear data
% save dataEllipse;