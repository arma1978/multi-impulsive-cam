clearvars; close all; clc; format longG;
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times');
set(0,'defaultUicontrolFontName','Times');
set(0,'defaultUitableFontName','Times');
set(0,'defaultTextFontName','Times');
set(0,'defaultUipanelFontName','Times');

projectRoot = fileparts(fileparts(mfilename('fullpath')));
cd(projectRoot);
addpath(fullfile(projectRoot,'src','routines'));
addpath(fullfile(projectRoot,'src','io'));
addpath(fullfile(projectRoot,'src','frames'));
addpath(fullfile(projectRoot,'utils'));
addpath(fullfile(getenv('HOME'),'Dropbox/Work/mosek/11.1/toolbox/r2022b'));

lowerPer = 5;
upperPer = 95;
%%
disc = 100;
sim = 1;
switch sim
    case 1
        tag = 'miss_distance';
        name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_J2miss_distance.mat');
    case 2
        tag = 'risk';
        name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_J2self_risk_Chan.mat');
    case 3
        tag = 'max_risk';
        name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_J2self_max_risk_estimate.mat'); 
end
load(name);

A = [out.majorIter];
[a,b] = find(A>2);

%% calculating the eccentricity
for ind = 1:length(out)
    out(ind).cVel = norm(out(ind).rrd(4:6)-out(ind).rrs(4:6));
end

%  A = [out.nImp];
%  [a,b] = find(A>165);

%[a,b] = find([out.cVel]<0.100);
[a,b] = find([out.ind]==928);

load(fullfile(projectRoot,'input','data.mat'))
B = data.Variables;
for ind = 1:size(data,1)
    rrd(1:3,ind) = B(ind,15:17)';
    rrd(4:6,ind) = B(ind,18:20)';
    rrs(1:3,ind) = B(ind,3:5);
    rrs(4:6,ind) = B(ind,6:8);
end

for ind = 1:length(b)
    [aa,bb] = min(abs(rrd(1,:)-out(b(ind)).rrd(1)));
    iddata(ind) = bb;
end

%% low relative velocity is 644
% 644   865   751   633   519   591   805   746
