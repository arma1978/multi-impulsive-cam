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
param.disc = 1111;
param.dvScl  = 1e6; %scaling for NLP
param.objScl = 1e2; % scaling for NLP
mass = 300; % kg
Tmax = 30e-3; %1000e-3; %15e-3; % 50e-3; %1000e-3;%13.0007e-3; %millinew

%% user input data
param.prop = 'J2';  % J2
param.nOrbits = 8; % time to maneuver in terms of period of sc
param.ind = 10; %1112; 1112          53         889         685         958        1026
param.order = 2;
param.nvar  = 7;
param.dt = 60; % in seconds
param.maxImpNum = 2000; % maximum number of impulses
param.dvMax = Tmax /mass*param.dt/1000;     % single dv km/s
PcLim = 1e-6;
maxPcLim = 1e-4;
minDistLim = 2;
sim = 3;

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
param.outfile = fullfile(outputDir,'ResultsNLP','sol');
param.fname   = fullfile(runtimeDir,'maps.dat');
param.projectRoot = projectRoot;
param.runtimeDir = runtimeDir;
if ~exist(fileparts(param.outfile),'dir'), mkdir(fileparts(param.outfile)); end

%% extract the test case and prepare the file for multi-impulse.
param.data = data; 
param = selectTestCaseFiltered(param);

%% run the main code (system for windows)
if strcmp(param.prop,'Kep')
    status = system(fullfile(projectRoot,'bin','runPropMultiMapsKep'));
    assert(status==0,'runPropMultiMapsKep execution failed.');
    dum = load(fullfile(runtimeDir,'xs0Kep.dat'));
    param.xx0s = dum(1:6);
else
    status = system(fullfile(projectRoot,'bin','runPropMultiMapsFullP'));
    assert(status==0,'runPropMultiMapsFullP execution failed.');
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
    Pb_2D0, maxPc0] = validate(zeros(nOptVarDvComp,1),DAx,param);
xx0 = [param.xx0s xx0; dvv0 zeros(3,1)];

%% convex solution
[out,xxOptLin,param] = linConvexSolve(param);

%% postporcess of the mosek solution
[sqrMahalanobisLin,tcaLin,dvLin,dvvLin,xxLin,...
    xxNoManLin,xxsEncLin,xxdEncLin,rrb2dLin,...
    Pb_2DLin, maxPcLin] = validate(xxOptLin,DAx,param);

toRTNLin = RTN(xxLin);
for ind = 1:size(dvvLin,2)
    dvvRTNLin(:,ind) = toRTNLin(:,:,ind)*dvvLin(:,ind);
    coeLin(:,ind) = pv2po(xxLin(1:3,ind), xxLin(4:6,ind), param.mu);
end

%% all nonlinear minimisation of dv, maha as constraint
objfnc  = @(xx) obj(xx, param);
nncofnc = @(xx) nncon(xx, DAx, param);

options = optimoptions('fmincon', 'Display', 'iter',...%'final-detailed',...
    'FiniteDifferenceStepSize', sqrt(eps), 'FiniteDifferenceType', 'central',...
    'Algorithm', 'sqp', 'MaxIterations', 10000, 'MaxFunctionEvaluations', 100000,...
    'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient',true,...
    'CheckGradients', false, 'UseParallel', true, 'ScaleProblem',true);

uub   =   param.dvMax*ones(param.ndv*3,1)*param.dvScl;
llb   =  -param.dvMax*ones(param.ndv*3,1)*param.dvScl;
tic
%xxOptLin = 1e-14*rand(size(xxOptLin));
xxOpt = fmincon(objfnc, xxOptLin*param.dvScl, [],[], [],[],llb, uub, nncofnc, options);
toc
xxOpt = xxOpt/param.dvScl;

%% postprocessing of fmincon 
[sqrMahalanobis,tca,dv,dvv,xx,xxNoMan,xxsEnc,xxdEnc,rrb2d, Pb_2D, maxPc] = validate(xxOpt,DAx,param);
%toRTN = RTN(xx);
toRTN = frenet(xx);

for ind = 1:size(dvv,2)
    dvvRTN(:,ind) = toRTN(:,:,ind)*dvv(:,ind);
    coe(:,ind) = pv2po(xx(1:3,ind), xx(4:6,ind), param.mu);
end
T = 2*pi*sqrt(coe(1,1:end-1).^3/param.mu);
tt = param.nOrbits-cumsum(param.dt./T);

%% writing summary 
fprintf('Convergence summary: \n');
fprintf('SOCP %i iterations  \n', size(out,2))
fprintf('SOCP delta v = %.4f mm/s \n', sum(dvLin)*km^2)
fprintf('SQP delta v = %.4f mm/s \n', sum(dv)*km^2)
fprintf('target squared Mahalanobis distance = %.4f \n', param.sqrMahaMin)
fprintf('SOCP squared Mahalanobis distance = %.4f \n', sqrMahalanobisLin)
fprintf('SQP Squared Mahalanobis distance = %.4f \n', sqrMahalanobis)
fprintf('initial tca = %.4f \n', param.t2TCA)
fprintf('SOCP  = %.4f \n', tcaLin)
fprintf('SQP tca   = %.4f \n', tca)

%
%% traj
figure()
quiver3(xx(1,:),xx(2,:),xx(3,:),dvv(1,:), dvv(2,:), dvv(3,:),'r', 'LineWidth', 2);
hold on
plot3(xxNoMan(1,1:end-1),xxNoMan(2,1:end-1),xxNoMan(3,1:end-1),'b', 'LineWidth', 2);
plot3(xxsEnc(1), xxsEnc(2), xxsEnc(3), 'rs', 'LineWidth', 2);
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal
grid on
box on

%% dv
figure(); hold on; grid on; box on
plot(tt,dvvRTN(1,1:end-1)*km, 'LineWidth', 2);
plot(tt,dvvRTN(2,1:end-1)*km, 'LineWidth', 2);
plot(tt,dvvRTN(3,1:end-1)*km, 'LineWidth', 2);
plot(tt,dv(1:end-1)*km, 'LineWidth', 2);
xlabel('Time to CA [rev]')
ylabel('$\Delta v$ [m/s]')
set(gca,'xdir','reverse')
legend('T', 'N', 'B', '\Delta v')

%% dv
figure(); hold on; grid on; box on
semilogy(tt,dv(1:end-1)*km,'LineWidth', 2)
xlabel('Time to CA [rev]')
ylabel('dv [m/s]')
set(gca,'xdir','reverse')

%% dv
figure(); hold on; grid on; box on
yyaxis right
plot(tt,dvvRTN(1,1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', [0, 0.4470, 0.7410], 'linestyle', '-');
plot(tt,dvvRTN(2,1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], 'linestyle', '-');
plot(tt,dvvRTN(3,1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', [0.9290, 0.6940, 0.1250], 'linestyle', '-');
plot(tt,dv(1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', 'k', 'linestyle', '-');
xlabel('Time to CA [rev]',  'color', 'k')
ylabel('Thrust [mN]',  'color', 'k')
set(gca,'xdir','reverse', 'ycolor', 'k')
legend('T', 'N', 'B', '\Delta v','location', 'northwest')
ylimmax = max([dvvRTN(1,1:end-1)*km/param.dt*mass*km ...
    dvvRTN(2,1:end-1)*km/param.dt*mass*km ...
    dvvRTN(3,1:end-1)*km/param.dt*mass*km ...
    dv(1:end-1)*km/param.dt*mass*km]);
ylimmin = min([dvvRTN(1,1:end-1)*km/param.dt*mass*km ...
    dvvRTN(2,1:end-1)*km/param.dt*mass*km ...
    dvvRTN(3,1:end-1)*km/param.dt*mass*km ...
    dv(1:end-1)*km/param.dt*mass*km]);

ylim(1.1*[ylimmin ylimmax])
yyaxis left
plot(tt,dvvRTN(1,1:end-1)*km^2, 'LineWidth', 2, 'color', [0, 0.4470, 0.7410], 'linestyle', '-');
plot(tt,dvvRTN(2,1:end-1)*km^2, 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], 'linestyle', '-');
plot(tt,dvvRTN(3,1:end-1)*km^2, 'LineWidth', 2, 'color', [0.9290, 0.6940, 0.1250], 'linestyle', '-');
plot(tt,dv(1:end-1)*km^2, 'LineWidth', 2, 'color', 'k', 'linestyle', '-');
xlabel('Time to CA [rev]', 'color', 'k')
ylabel('$\Delta v$ [mm/s]',  'color', 'k')
set(gca,'xdir','reverse', 'ycolor', 'k')
legend('T', 'N', 'B', '||\cdot||','location', 'northeast')
ylim(1.1*[ylimmin ylimmax]./mass*param.dt)

%% dv
figure(); hold on; grid on; box on
semilogy(tt,dv(1:end-1)*km,'LineWidth', 2)
xlabel('Time to CA [rev]')
ylabel('dv [m/s]')
set(gca,'xdir','reverse')

%% save workspace and output file
clear traindata_augmented_cov
fnamemat = strcat([param.outfile,param.prop,'_case_',num2str(param.ind),...
    '_time2man_', num2str(param.nOrbits*100),'_maxCollPro_', num2str(-log10(param.PcMax))]);
save(fnamemat);