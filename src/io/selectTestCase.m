function [param, covsRTN, covdRTN] = selectTestCase(param)
% selectTestCase  Extract one conjunction row from the legacy catalogue format.
%
% This variant is used by mainLargeSimSCVX/mainMILP on conjunctions.mat.
% It writes encounter states to disk, runs backward propagation, and writes
% inputFullP.dat for the forward multi-impulse propagator.

km = 1000;
if isfield(param,'projectRoot')
    projectRoot = param.projectRoot;
else
    projectRoot = pwd;
end
if isfield(param,'runtimeDir')
    runtimeDir = param.runtimeDir;
else
    runtimeDir = fullfile(projectRoot,'runtime');
end
B = param.data;
param.TestCase = B(param.ind,:);
B = B(:,108:end);
B = B.Variables; 
param.bodySize = (param.TestCase{1,'c_span'}+param.TestCase{1,'t_span'})/1000; % km

%% extract data for the spacecraft all in in m m/s converted into km km/s
pos = B(param.ind,1:6);
pos(1) = pos(1)/km;
Ts = 2*pi*sqrt(pos(1)^3/param.mu);
param.t2TCA = Ts*param.nOrbits;
pos(3:6) = pos(3:6)*pi/180;
[rrs, vvs] = po2pv(pos, param.mu);

covsRTN = [B(param.ind,17) B(param.ind,20)  B(param.ind,21);
          B(param.ind,20)  B(param.ind, 18) B(param.ind,22);
          B(param.ind,21)  B(param.ind, 22)  B(param.ind,19)]./km^2;
%% writting to file the data of the spacecraft in order to backpropagate 
% to simulate an encounter
fid = fopen(fullfile(runtimeDir,'xsTCA.dat'), 'w');
fprintf(fid, '%40.12f\n', rrs(1));
fprintf(fid, '%40.12f\n', rrs(2));
fprintf(fid, '%40.12f\n', rrs(3));
fprintf(fid, '%40.12f\n', vvs(1));
fprintf(fid, '%40.12f\n', vvs(2));
fprintf(fid, '%40.12f\n', vvs(3));
fprintf(fid, '%40.12f\n', param.t2TCA);
fprintf(fid, '%40.12f\n', param.dt);
fclose(fid);

%% getting the data of the debris, same as before
pod = B(param.ind,7:12);
pod(1) = pod(1)/km;
pod(3:6) = pod(3:6)*pi/180; 
[rrd, vvd] = po2pv(pod, param.mu);

covdRTN = [B(param.ind,23) B(param.ind,26) B(param.ind,27);
          B(param.ind,26)  B(param.ind,24) B(param.ind,28);
          B(param.ind,27)  B(param.ind,28) B(param.ind,25)]./km^2;
      
%% writting to file the data of the debris in order to backpropagate 
% to simulate an encounter

fid = fopen(fullfile(runtimeDir,'xdTCA.dat'), 'w');
fprintf(fid, '%40.12f\n', rrd(1));
fprintf(fid, '%40.12f\n', rrd(2));
fprintf(fid, '%40.12f\n', rrd(3));
fprintf(fid, '%40.12f\n', vvd(1));
fprintf(fid, '%40.12f\n', vvd(2));
fprintf(fid, '%40.12f\n', vvd(3));
fprintf(fid, '%40.12f\n', param.t2TCA);
fprintf(fid, '%40.12f\n', param.dt);
fclose(fid);

%% Prepare the input file for the forward multi-impulse collision avoidance
%!./runStateBackProp
status = system(fullfile(projectRoot,'bin','stateBackProp'));
assert(status==0,'stateBackProp execution failed.');

%param.xx0d = [rrd; vvd];
param.PdRTN = covdRTN;
%param.xx0s = [rrs; vvs];
param.PsRTN = covsRTN;
%param = rmfield(param,'data'); % remove data as too heavy
param.coed=pod';
param.coes=pos';
param.Ts = Ts;

%% write full input file 
fid = fopen(fullfile(runtimeDir,'inputFullP.dat'), 'w');
fprintf(fid, '%2i\n', param.order);
fprintf(fid, '%40.12f\n', param.dt);
fprintf(fid, '%40.12f\n', param.bodySize);
fprintf(fid, '%40.12f\n', param.PsRTN(1,1));
fprintf(fid, '%40.12f\n', param.PsRTN(2,2));
fprintf(fid, '%40.12f\n', param.PsRTN(3,3));
fprintf(fid, '%40.12f\n', param.PsRTN(1,2));
fprintf(fid, '%40.12f\n', param.PsRTN(1,3));
fprintf(fid, '%40.12f\n', param.PsRTN(2,3));
fprintf(fid, '%40.12f\n', param.PdRTN(1,1));
fprintf(fid, '%40.12f\n', param.PdRTN(2,2));
fprintf(fid, '%40.12f\n', param.PdRTN(3,3));
fprintf(fid, '%40.12f\n', param.PdRTN(1,2));
fprintf(fid, '%40.12f\n', param.PdRTN(1,3));
fprintf(fid, '%40.12f\n', param.PdRTN(2,3));
fprintf(fid, '%4i\n', param.maxImpNum);
fclose(fid);
