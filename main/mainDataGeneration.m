% =========================================================================
% mainDataGeneration.m
% Build derived conjunction datasets from conjunctions.mat for campaign runs.
% Requires local input/conjunctions.mat (not versioned in this repository).
% =========================================================================

clearvars; 
close all; 
clc

projectRoot = fileparts(fileparts(mfilename('fullpath')));
cd(projectRoot);

load(fullfile(projectRoot,'input','conjunctions.mat'));
PcLim = 1e-6;
maxPcLim = 1e-4;
minDistLim = 2;
km = 1000;

dataRev = traindata_augmented_cov;

cond1 = dataRev.self_det_C_B<1e16;
cond2 = dataRev.self_det_C_B>1e5;
cond3 = dataRev.miss_distance<minDistLim*km;
cond4 = dataRev.self_max_risk_estimate>log10(maxPcLim);
cond5 = dataRev.self_risk_Chan>log10(PcLim);
dataRev = dataRev(cond1&cond2&cond3&cond4&cond5,:);

save(fullfile(projectRoot,'input','dataRev.mat'),'dataRev')