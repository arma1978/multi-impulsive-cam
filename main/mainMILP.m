% =========================================================================
% mainMILP.m
% Solve the multi-impulse collision-avoidance problem using the MILP model.
% Requires local input/conjunctions.mat (not versioned in this repository).
% =========================================================================

clearvars; close all; clc; format longG;
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
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
load(fullfile(projectRoot,'input','conjunctions.mat'));
km = 1000;
param.mu = 398600.4418;
param.rE = 6378.137;
param.J2 = 1.08262668e-3;
param.Lsc = param.rE;
param.Vsc = sqrt(param.mu/param.Lsc);
param.Tsc = param.Lsc/param.Vsc;
param.disc = 301;

%% user input data
param.prop = 'J2'; %J2
param.nOrbits = 12; % time to maneuver in terms of period of sc
param.ind = 1; % which test case
param.order = 4;
param.nvar  = 7;
param.dt = 60; % in seconds
param.maxImpNum = 400; % maximum number of impulses
param.dvMax = 0.06/1000*param.dt/1000;     % single dv km/s
param.sqrMahaMin = []; % if this is not specified is calculated from PcMax
param.PcMax = 1e-6;
param.Pc      = [];   % target Pc ([] means use PcMax instead)
param.minDist = [];   % minimum miss distance [km] ([] means use PcMax/Pc)

%% parameters setting
param.outfile = fullfile(outputDir,'ResultsMILP','sol');
param.fname   = fullfile(runtimeDir,'maps.dat');
param.projectRoot = projectRoot;
param.runtimeDir = runtimeDir;
if ~exist(fileparts(param.outfile),'dir'), mkdir(fileparts(param.outfile)); end

%% extract the test case and prepare the file for multi-impulse.
param.data = traindata_augmented_cov;
param.data = sortrows(param.data,'max_risk_estimate', 'descend');
param = selectTestCase(param);

%% run the main code (system for windows)
if strcmp(param.prop,'Kep')
    status = system(fullfile(projectRoot,'bin','statePropMultiMapsKep'));
    assert(status==0,'statePropMultiMapsKep execution failed.');
    dum = load(fullfile(runtimeDir,'xs0Kep.dat'));
    param.xx0s = dum(1:6);
else
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
nOptVar = N*6 + param.disc -1;
param.nOptVar = nOptVar;
param.ndv = (nOptVar-param.disc+1)/6;

%% pre-calculating relevant quantities
param = refTraj(DAx,param);
param = linearMaps(DAx,param);
param = ellipseTargets(param);

[sqrMaha0,tca0,dv0,dvv0,xx0,...
    xxNoMan0,xxsEnc0,xxdEnc0,rrb2d0,...
    Pb_2D0, maxPc0] = validate(zeros(param.ndv*3,1),DAx,param);

%% MILP solver
[f,A,b,Aeq,beq,intcon,uubLin,llbLin,param] = linProgProblem(param);
solLin = intlinprog(f,intcon,A,b,[],[],llbLin,uubLin);
xxOptLin = solLin(1:param.ndv*3) - solLin(param.ndv*3+1:param.ndv*6);

%% postprocessing  of the solution 
[sqrMahalanobisLin, tcaLin, dvLin, dvvLin, ...
    rrvvLin, rrvvNoManLin, rrvvsEncLin, rrvvdEncLin, ...
    rrb2dLin, Pb_2DLin, maxPcLin] = validate(xxOptLin,DAx,param);

toRTN = RTN(rrvvLin);
for ind = 1:size(dvvLin,2)
    dvvRTN(:,ind) = toRTN(:,:,ind)*dvvLin(:,ind);
    coe(:,ind) = pv2po(rrvvLin(1:3,ind), rrvvLin(4:6,ind), param.mu);
end
T = 2*pi*sqrt(coe(1,1:end-1).^3/param.mu);
tt = param.nOrbits-cumsum(param.dt./T);

%% output summary
fprintf('Mosek convergence summary: \n');
fprintf('delta v 1st round = %.4f m/s \n', sum(dvLin)*1000)
fprintf('squared Mahalanobis distance 1st round = %.4f \n', sqrMahalanobisLin)
fprintf('target squared Mahalanobis distance = %.4f \n', param.sqrMahaMin)
fprintf('initial tca = %.4f \n', param.t2TCA)
fprintf('tca 1st round  = %.4f \n', tcaLin)

%% plotting trj
figure();
quiver3(rrvvLin(1,:),rrvvLin(2,:),rrvvLin(3,:), dvvLin(1,:), dvvLin(2,:), dvvLin(3,:),'r', 'LineWidth', 2);
hold on; 
plot3(rrvvNoManLin(1,1:end-1),rrvvNoManLin(2,1:end-1),rrvvNoManLin(3,1:end-1),'b','LineWidth', 2);
plot3(rrvvsEncLin(1), rrvvsEncLin(2), rrvvsEncLin(3), 'rs', 'LineWidth', 2);
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal
grid on; 
box on

%% dv
figure(); hold on; grid on; box on
plot(tt,dvvRTN(1,1:end-1)*km, 'LineWidth', 2);
plot(tt,dvvRTN(2,1:end-1)*km, 'LineWidth', 2);
plot(tt,dvvRTN(3,1:end-1)*km, 'LineWidth', 2);
plot(tt,dvLin(1:end-1)*km, 'LineWidth', 2);
xlabel('Time to conjuction [rev]')
ylabel('$\Delta v$ [m/s]')
set(gca,'xdir','reverse')
legend('R', 'T', 'N', 'mag')

% Bplane
figure()
plot(param.xx2Target(2,:), param.xx2Target(1,:), 'LineWidth', 2)
hold on 
plot(rrb2dLin(2,:),rrb2dLin(1,:), 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k')
%plot(param.recTargetPoint(2,:),param.recTargetPoint(1,:),'g','LineWidth', 2)
ylabel('$\xi$ [km]')
xlabel('$\zeta$ [km]')
grid on 
box on 
axis equal

%% save workspace and output file
clear traindata_augmented_cov
fnamemat = strcat([param.outfile,param.prop,'_case_',num2str(param.ind),...
    '_time2man_', num2str(param.nOrbits*100),'_maxCollPro_', num2str(-log10(param.PcMax))]);
save(fnamemat);