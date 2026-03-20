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

for const = 1:3
    if const==1
        tag = 'max_risk';
        name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_J2self_max_risk_estimate.mat');
    elseif const==2
        tag = '$P_{C}$ [-]';
        name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_J2self_risk_Chan.mat');
    else
        tag = 'miss_distance';
        name = fullfile(projectRoot,'output','ResultsLargeSimSCVX','Rev_J2miss_distance.mat');
    end
    load(name);
    BsqrMaha(:,const) = [out.sqrMahalanobisRef]';
    Bdv(:,const) = [out.dv]';
    BnImp(:,const) = [out.nImp]';
    BminrelDist(:,const) = [out.minRelDistRef]';
    BmaxPc(:,const) = [out.maxPcRef]';
    BPc(:,const) = [out.PcRef]';
end

 
%% plotting the results
res = '-r500';
type = '-dpng';
path = fullfile(projectRoot,'output','Figures');
if ~exist(path,'dir'), mkdir(path); end

figure
nameFig = 'BoxPlot_dv';
boxplot(Bdv*1000, 'labels',{'$P_{C,\max}$ [-]', '$P_{C}$ [-]', '$\Delta r_{CA}$'}, 'whisker', 1.5, 'symbol', '' )
ylabel('$\Delta v$ [m/s]')
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
savePlot(nameFig, path, type, res);

figure
nameFig = 'BoxPlot_nImp';
boxplot(BnImp, 'labels',{'$P_{C,\max}$ [-]', '$P_{C}$ [-]', '$\Delta r_{CA}$'}, 'whisker', 1.5, 'symbol', '' )
ylabel('Impulse $\#$ [-]')
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
savePlot(nameFig, path, type, res);

figure
nameFig = 'BoxPlot_relDist';
boxplot(BminrelDist, 'labels',{'$P_{C,\max}$ [-]', '$P_{C}$ [-]', '$\Delta r_{CA}$'},  'whisker', 1.5, 'symbol', '' )
ylabel('$\Delta r_{CA}$ [km]')
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
savePlot(nameFig, path, type, res);

figure
nameFig = 'BoxPlot_Pcmax';
boxplot(BmaxPc, 'labels',{'$P_{C,\max}$ [-]', '$P_{C}$ [-]', '$\Delta r_{CA}$'},  'whisker', 1.5, 'symbol', '' )
ylabel('$P_{C,\max}$ [-]')
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
savePlot(nameFig, path, type, res);

figure
nameFig = 'BoxPlot_Pc';
boxplot(BPc, 'labels',{'$P_{C,\max}$ [-]', '$P_{C}$ [-]', '$\Delta r_{CA}$'},  'whisker', 1.5, 'symbol', '' )
ylabel('$P_{C}$ [-]')
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
savePlot(nameFig, path, type, res);
