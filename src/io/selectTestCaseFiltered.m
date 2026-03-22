function [param, covsRTN, covdRTN] = selectTestCaseFiltered(param)
% selectTestCaseFiltered  Extract and prepare a single conjunction test case.
%
%   Reads the conjunction catalogue stored in param.data, extracts the row
%   at index param.ind, computes the manoeuvre window duration, writes the
%   spacecraft and debris initial conditions to disk for the C++ propagator,
%   runs the backward propagator (stateBackProp) to obtain the initial state
%   at the start of the manoeuvre window, and writes the full propagation
%   input file (inputFullP.dat).
%
% Inputs: param  - struct with fields: data, ind, prop, nOrbits, dt,
%                   maxImpNum, order, bodySize, mu, J2, rE
% Outputs: param    - updated struct (adds: t2TCA, PsRTN, PdRTN, coes, coed, Ts)
%          covsRTN  - 3x3 chaser covariance in RTN frame  [km^2]
%          covdRTN  - 3x3 debris  covariance in RTN frame  [km^2]

B = param.data;
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
param.TestCase = B(param.ind, :);
B = B.Variables;
param.bodySize = B(param.ind, 2);  % combined hard-body radius [km]

%% Extract spacecraft state at TCA  [km, km/s]
rrs = B(param.ind, 3:5);
vvs = B(param.ind, 6:8);

%% Compute orbital period and manoeuvre window length
if strcmp(param.prop, 'J2')
    % Use J2-corrected nodal period
    E  = pv2po(rrs', vvs', param.mu);
    Ts = periods(E(1), E(2), E(3), E(5), param.J2, param.rE, param.mu);
else
    % Two-body Keplerian period
    en = 0.5*norm(vvs)^2 - param.mu/norm(rrs);
    a  = -param.mu / (2*en);
    Ts = 2*pi * sqrt(a^3 / param.mu);
end

param.t2TCA = Ts * param.nOrbits;  % total manoeuvre window [s]

covsRTN = [B(param.ind,9) B(param.ind,12)  B(param.ind,13);
    B(param.ind,12)  B(param.ind, 10) B(param.ind,14);
    B(param.ind,13)  B(param.ind, 14)  B(param.ind,11)];
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
rrd = B(param.ind,15:17);
vvd = B(param.ind,18:20);

covdRTN = [B(param.ind,21) B(param.ind,24) B(param.ind,25);
    B(param.ind,24)  B(param.ind,22) B(param.ind,26);
    B(param.ind,25)  B(param.ind,26) B(param.ind,23)];

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

%% this preparaes the input file for the forward multi-impulse collision avoidance
%!./runStateBackProp
status = system(fullfile(projectRoot,'bin','stateBackProp'));
assert(status==0,'stateBackProp execution failed.');

%param.xx0d = [rrd; vvd];
param.PdRTN = covdRTN;
%param.xx0s = [rrs; vvs];
param.PsRTN = covsRTN;
%param = rmfield(param,'data'); % remove data as too heavy
param.coed=pv2po(rrd',vvd',param.mu);
param.coes=pv2po(rrs',vvs',param.mu);
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
