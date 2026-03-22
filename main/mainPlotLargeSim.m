% =========================================================================
% mainPlotLargeSim.m
% Generate aggregate plots for large-simulation campaign results.
% Requires generated MAT files in output/ResultsLargeSimSCVX/.
% =========================================================================

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
sim = 2;
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

A = [out.sqrMahalanobisLin];
B = [out.sqrMahalanobisRef];
close all
plot(sort(abs(B-A)./B*100))

%% plotting the results
res = '-r500';
type = '-dpng';
path = fullfile(projectRoot,'output','Figures');
if ~exist(path,'dir'), mkdir(path); end

%% 1
data= sortrows(data,'Pc_approx', 'descend');
[a,b] = sort([out.Pc0],'descend');
out = out(b);

figure(); hold on; grid on; box on
name = 'PcData';
nameFig = strcat(tag,name);
plot(data.Pc,'LineWidth', 2);
set(gca, 'yScale', 'log')
xlabel('Test case \# [-]');
ylabel('$P^*_C$ [-]');
savePlot(nameFig, path, type, res);

%%2
data= sortrows(data,'Pc_max', 'descend');
[a,b] = sort([out.maxPc0],'descend');
out = out(b);

figure(); hold on; grid on; box on
name = 'PcmaxData';
nameFig = strcat(tag,name);
plot([1:size(data,1)],data.Pc_max,'LineWidth', 2);
set(gca, 'yScale', 'log')
xlabel('Test case \# [-]');
ylabel('$P^*_{C,\max}$ [-]');
savePlot(nameFig, path, type, res);

%% 3
data= sortrows(data,'d^* [km]', 'ascend');
[a,b] = sort([out.minRelDist0],'descend');
out = out(b);

figure(); hold on; grid on; box on
name = 'minDist';
nameFig = strcat(tag,name);
semilogy([1:size(data,1)],data.("d^* [km]"),'LineWidth', 2);
xlabel('Test case \# [-]');
ylabel('$\Delta r_{CA}^*$ [km]');
savePlot(nameFig, path, type, res);

%% plotting
name = 'sqrMaha';
nameFig = strcat(tag,name);
figure(); hold on; grid on; box on
histogram(rmoutliers([out.sqrMahalanobisRef], 'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$d_{CA}^2$ [km$^2$]')
ylabel('Probability')
savePlot(nameFig, path, type, res);

%% plotting
name = 'dv';
nameFig = strcat(tag,name);
figure(); hold on; grid on; box on
histogram(rmoutliers([out.dv], 'percentile', [lowerPer upperPer])*km^2,disc, 'Normalization', 'probability')
xlabel('$\Delta v$ [mm/s]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
name = 'nImp';
nameFig = strcat(tag,name);
figure(); hold on; grid on; box on
histogram(rmoutliers([out.nImp],'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('Number of impulses [-]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
name = 'dtca';
nameFig = strcat(tag,name);
figure(); hold on; grid on; box on
histogram(rmoutliers([out.tcaRef]-[out.tca0], 'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$\Delta t_{CA}$ [s]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'minRelDist';
nameFig = strcat(tag,name);
histogram(rmoutliers([out.minRelDistRef],'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$\Delta r_{CA}$ [km]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'maxPc';
nameFig = strcat(tag,name);
histogram(rmoutliers([out.maxPcRef],'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$P_{C,\max}$ [-]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'Pc';
nameFig = strcat(tag,name);
histogram(rmoutliers([out.PcRef],'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$P_{C}$ [-]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'compTime';
nameFig = strcat(tag,name);
histogram(rmoutliers([out.compTime]+([out.majorIter]+1)*0.05,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('Computational Time [s]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'ErrorOnDistance';
nameFig = strcat(tag,name);
histogram(rmoutliers(abs([out.minRelDistLin] - [out.minRelDistRef])*km,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('Error on $\Delta r_{CA}$ [m]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);


%% plotting
figure(); hold on; grid on; box on
name = 'ErrorOnDistance';
nameFig = strcat(tag,name);
histogram(rmoutliers(abs([out.minRelDistLin] - [out.minRelDistRef])*km,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('Error on $\Delta r_{CA}$ [m]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'ErrorOnPcMax';
nameFig = strcat(tag,name);
histogram(rmoutliers(abs([out.maxPcRef] - [out.maxPcLin]),'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('Error on $P_{C,\max}$ [-]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'ErrorOnPc';
nameFig = strcat(tag,name);
histogram(rmoutliers(abs([out.PcRef] - [out.PcLin]),'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('Error on $P_{C}$ [-]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);

%% plotting
figure(); hold on; grid on; box on
name = 'ErrorLin';
nameFig = strcat(tag,name);
if isfield(out,'dvLin')
    errLin = abs([out.dvLin] - [out.dv]) * km^2;
else
    % Fallback for result sets that do not store linear-round delta-v.
    errLin = abs([out.dv]) * km^2;
end
histogram(rmoutliers(errLin,'percentile', [lowerPer upperPer]),disc, 'Normalization', 'probability')
xlabel('$||\Delta {\bf v} - \Delta {\bf v}_{Lin}||$ [mm/s]')
ylabel('Probability [-]')
savePlot(nameFig, path, type, res);


%% calculating the eccentricity
for ind = 1:length(out)
    podumd = pv2po(out(ind).rrd(1:3), out(ind).rrd(4:6),3.98e5);
    podums = pv2po(out(ind).rrs(1:3), out(ind).rrs(4:6),3.98e5);
    out(ind).eccd = podumd(2);
    out(ind).eccs = podums(2);
end


%data= sortrows(data,'Pc', 'descend');
%data{:,'ID'}= [1:size(data.ID,1)]';
%data.('tca [s]') = [];
%writetable(data,'data.dat')
%save data.mat data