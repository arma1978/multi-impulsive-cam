% =========================================================================
% mainLargeSimSCVX.m
% Run the SCVX pipeline over a large conjunction catalogue and save results.
% Requires local input/conjunctions.mat (not versioned in this repository).
% =========================================================================

clearvars -except quickTestMaxCases; close all; clc; format longG;
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

%% input
PcLim = 1e-6;
maxPcLim = 1e-4;
minDistLim = 2;
km = 1000;

%% DATA %% all units are km, s,
load(fullfile(projectRoot,'input','conjunctions.mat'));
param.mu = 398600.4418;
param.rE = 6378.137;
param.J2 = 1.08262668e-3;
param.Lsc = param.rE;
param.Vsc = sqrt(param.mu/param.Lsc);
param.Tsc = param.Lsc/param.Vsc;
param.disc = 1111;

%% user input data
param.prop = 'Kep';  % J2
param.nOrbits = 2; % time to maneuver in terms of period of sc
param.order = 2;
param.nvar  = 7;
param.dt = 60; % in seconds
param.maxImpNum = 170; % maximum number of impulses
param.dvMax = 1e-4*param.dt/1000;     % single dv km/s

%% parameters setting
param.outfile = fullfile(outputDir,'ResultsLargeSimSCVX','Rev_');
param.fname   = fullfile(runtimeDir,'maps.dat');
param.projectRoot = projectRoot;
param.runtimeDir = runtimeDir;
if ~exist(fileparts(param.outfile),'dir'), mkdir(fileparts(param.outfile)); end

sim = 1;
switch sim
    case 1
        obj = 'miss_distance';
        order = 'ascend';
        param.minDist = minDistLim; % rivedere!! provide this for min distance
        param.PcMax = maxPcLim;
        param.Pc = [];
    case 2
        obj = 'self_risk_Chan';
        order = 'descend';
        param.Pc = PcLim;
        param.minDist = [];
        param.PcMax = [];
    case 3
        obj = 'self_max_risk_estimate';
        order = 'descend';
        param.PcMax = maxPcLim;
        param.minDist = [];
        param.Pc = [];
end

%% extract the test case and prepare the file for multi-impulse.
param.data = traindata_augmented_cov;
param.data = sortrows(param.data,obj, order);

%filtering data
cond1 = param.data.self_det_C_B<1e16;
cond2 = param.data.self_det_C_B>1e5;
cond3 = param.data.miss_distance<minDistLim*km;
cond4 = param.data.self_max_risk_estimate>log10(maxPcLim);
cond5 = param.data.self_risk_Chan>log10(PcLim);
param.data = param.data(cond1&cond2&cond3&cond4&cond5,:);

% Set quickTestMaxCases (e.g., 3) before running for a short smoke test.
if ~exist('quickTestMaxCases','var') || isempty(quickTestMaxCases)
    quickTestMaxCases = size(param.data,1);
end
quickTestMaxCases = min(quickTestMaxCases,size(param.data,1));
fprintf('Running %i/%i filtered test cases.\n',quickTestMaxCases,size(param.data,1));

figure(); hold on; grid on; box on
plot(param.data{:,obj},'LineWidth', 2);
xlabel('Test Case');
ylabel(obj);
line = [];

data  = array2table(zeros(0,33), ...
    'VariableNames', {'ID', 'tca [s]', 'R [km]', 'p_j2k_x [km]', 'p_j2k_y [km]', 'p_j2k_z [km]', 'p_j2k_vx [km/s]', 'p_j2k_vy [km/s]', 'p_j2k_vz [km/s]',...
    'p_c_rr  [km^2]', 'p_c_tt  [km^2]', 'p_c_nn  [km^2]', 'p_c_rt  [km^2]', 'p_c_rn  [km^2]', 'p_c_tn  [km^2]'...
    's_j2k_x [km]',   's_j2k_y [km]', 's_j2k_z [km]', 's_j2k_vx [km/s]', 's_j2k_vy [km/s]', 's_j2k_vz [km/s]',...
    's_c_rr  [km^2]', 's_c_tt  [km^2]', 's_c_nn  [km^2]', 's_c_rt  [km^2]', 's_c_rn  [km^2]', 's_c_tn  [km^2]',...
    'Pc', 'Pc_approx', 'Pc_max', 'd^* [km]', 'v^* [km/s]', 'd_m^2 [km^2]'});
figure()
for ind = 1:quickTestMaxCases
    % waitbar(ind/size(param.data,1));
    param.ind = ind;
    fprintf('Test case number %i/%i \n', ind,size(param.data,1));
    [param, covsRTN, covdRTN] = selectTestCase(param);
    
    %% run the main code (system for windows)
    if strcmp(param.prop,'Kep')
        %!./runPropMultiMapsKep
        status = system(fullfile(projectRoot,'bin','statePropMultiMapsKep'));
        assert(status==0,'statePropMultiMapsKep execution failed.');
        dum = load(fullfile(runtimeDir,'xs0Kep.dat'));
        param.xx0s = dum(1:6);
    else
        %!./runPropMultiMapsFullP
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
    tic;
    [outLin,xxOptLin,param] = linConvexSolve(param);
    out(ind).compTime = toc;
    [sqrMahalanobisLin, tcaLin, dvLin, dvvLin, ...
        rrvvLin, rrvvNoManLin, rrvvsEncLin, rrvvdEncLin, ...
        rrb2dLin, Pb_2DLin, maxPcLin, PcLin] = validate(xxOptLin,DAx,param);
    
    %% major iterations
    tolM = Inf;
    xx(:,1) = xxOptLin;
    tca(1) = tcaLin;
    pp(1) = param;
    jj=0;
    minorIter = size(outLin,2);
    outRef = struct('compTime',0);
    while (tolM>1e-6)
        jj = jj+1;
        [outRef,xx(:,jj+1),DAxRef,pp(jj+1)] = refConvex(xx(:,jj), tca(jj), pp(jj));
        out(ind).compTime = out(ind).compTime+sum([outRef.compTime]);
        [~, tca(jj+1)] = validate(xx(:,jj+1)-xx(:,jj),DAxRef,pp(jj+1));
        tolM = max(abs(xx(:,jj+1)-xx(:,jj)));
        minorIter = minorIter + size(outRef,2);
    end
    paramRef = pp(end);
    
    %% postprocessing
    [sqrMahalanobisRef, tcaRef, dvRef, dvvRef, ...
        rrvvRef, rrvvNoManRef, rrvvsEncRef, rrvvdEncRef, ...
        rrb2dRef, Pb_2DRef, maxPcRef, PcRef] = validate(xx(:,jj+1)-xx(:,jj),DAxRef,paramRef);
    
    dvv = [reshape(xx(:,jj+1),3,param.ndv) zeros(3,1)];
    dv = sqrt(sum(dvv.*dvv,1));
    
    %% output summary
    %% output summary
    fprintf('Mosek convergence summary: \n');
    fprintf('%i iterations for 1-st round, %i iterations for %i-th round  \n', size(outLin,2),size(outRef,2),jj+1)
    fprintf('delta v 1-st round = %.4f m/s \n', sum(dvLin)*1000)
    fprintf('delta v %i-th round = %.4f m/s \n', jj+1, sum(dv)*1000)
    fprintf('squared Mahalanobis distance 1-st round = %.4f \n', sqrMahalanobisLin)
    fprintf('squared Mahalanobis distance %i-th round = %.4f \n', jj+1, sqrMahalanobisRef)
    fprintf('target squared Mahalanobis distance = %.4f \n', param.sqrMahaMin)
    fprintf('initial tca = %.4f \n', param.t2TCA)
    fprintf('tca 1-st round  = %.4f \n', tcaLin)
    fprintf('tca %i-th round = %.4f \n',  jj+1, tcaRef)
    
    out(ind).ind = ind;
    out(ind).sqrMaha0 = sqrMaha0;
    out(ind).tca0 = tca0;
    out(ind).covsRTN = covsRTN;
    out(ind).covdRTN = covdRTN;
    out(ind).rrs = xxsEnc0;
    out(ind).rrd = xxdEnc0;
    out(ind).maxPc0 = maxPc0;
    out(ind).Pc0 = Pc0;
    out(ind).minRelDist0 = norm(rrb2d0);
    
    out(ind).sqrMahalanobisLin = sqrMahalanobisLin;
    out(ind).tcaLin = tcaLin;
    out(ind).maxPcLin = maxPcLin;
    out(ind).PcLin = PcLin;
    out(ind).dvLin = sum(dvLin);
    out(ind).minRelDistLin = norm(rrb2dLin);
    out(ind).rrb2dLin = rrb2dLin;
    
    
    out(ind).sqrMahalanobisRef = sqrMahalanobisRef;
    out(ind).tcaRef = tcaRef;
    out(ind).maxPcRef = maxPcRef;
    out(ind).PcRef = PcRef;
    out(ind).minRelDistRef = norm(rrb2dRef);
    out(ind).rrb2dRef = rrb2dRef;
    
    out(ind).iterFirst = size(outLin,2);
    out(ind).iterSecond = size(outRef,2);
    out(ind).dvv = dvv;
    out(ind).dv = sum(dv);
    out(ind).nImp = sum(dv>1e-2*param.dvMax);
    out(ind).PcAlfano = param.PcAlfano;
    out(ind).majorIter = jj;
    out(ind).minorIter = minorIter;
    out(ind).compTime = out(ind).compTime + sum([outRef.compTime]);
    line = [ind tca0 param.bodySize xxsEnc0' covsRTN(1,1) covsRTN(2,2) covsRTN(3,3) covsRTN(1,2) covsRTN(1,3) covsRTN(2,3)...
        xxdEnc0' covdRTN(1,1) covdRTN(2,2) covdRTN(3,3) covdRTN(1,2) covdRTN(1,3) covdRTN(2,3)...
        param.PcAlfano Pc0 maxPc0 norm(rrb2d0) norm(xxsEnc0(4:6)-xxdEnc0(4:6)) sqrMaha0];
    data{ind,:} = line;
    clear xx pp tca; 
end

%% filtering out data
out(([out.maxPc0]<maxPcLim)|([out.maxPc0]>1)) = [];
out([out.minRelDist0]>minDistLim) = [];
out([out.Pc0]<PcLim) = [];

cond1 = data.Pc_max<1;
cond2 = data.("d^* [km]")<minDistLim;
cond3 = data.Pc_max>maxPcLim;
cond4 = data.Pc>PcLim;

data = data(cond1&cond2&cond3&cond4,:);
data = sortrows(data,'d^* [km]', 'ascend');
data{:,'ID'}= [1:size(data.ID,1)]';
[a,b] = sort([out.minRelDist0],'ascend');
out = out(b);

%% save workspace and output file
clear traindata_augmentes_cov
fnamemat = strcat([param.outfile,param.prop,obj]);
save(fnamemat);