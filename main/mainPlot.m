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

name = fullfile(projectRoot,'output','ResultsSCVX','solJ2_ID_1_obj_max_risk_time2man_2_maxnImp_180.mat');
tag  = 'max_risk';
load(name);

%% plotting the results
res = '-r500';
type = '-dpng';
path = fullfile(projectRoot,'output','Figures');
if ~exist(path,'dir'), mkdir(path); end

%% plotting
figure(); hold on; grid on; box on
name = 'b-plane';
nameFig = strcat(tag,name);
plot(param.xx2Target(2,:),param.xx2Target(1,:),'k', 'LineWidth', 2)
plot(paramRef.xx2Target(2,:),paramRef.xx2Target(1,:),'b', 'LineWidth', 2)
dumT = [outLin.xx2Opt outRef.xx2Opt];
dumv = [outLin.dvTot outRef.dvTot]*1000;
for ind = 1:length(dumv)
    textstr{ind} = num2str(ind);
end
%text(dumT(2,:),dumT(1,:),textstr, 'FontSize', 20, 'FontName', 'Times')
scatter(dumT(2,:), dumT(1,:), 80, dumv','filled','MarkerEdgeColor', 'k')
dumT = [outLin.xxOnEllipse outRef.xxOnEllipse];
%plot(dumT(2,:), dumT(1,:), 'o', 'MarkerSize', 7, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5)
%text(dumT(2,:),dumT(1,:),textstr, 'FontSize', 20, 'FontName', 'Times')
%plot(param.rrb2Nom(2), param.rrb2Nom(1), 'ko', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
%plot(rrb2dRef(2), rrb2dRef(1), 'ko', 'MarkerFaceColor','g', 'MarkerEdgeColor', 'g')
ylabel('$\xi$ [km]')
xlabel('$\zeta$ [km]')
colorbar
title('$\Delta v$ [m/s]')
axis equal
savePlot(nameFig, path, type, res);

%% traj
figure()
name = 'trj';
nameFig = strcat(tag,name);
quiver3(rrvvRef(1,:),rrvvRef(2,:),rrvvRef(3,:),dvv(1,:), dvv(2,:), dvv(3,:),'r', 'LineWidth', 2);
hold on
plot3(rrvvNoManRef(1,1:end-1),rrvvNoManRef(2,1:end-1),rrvvNoManRef(3,1:end-1),'b', 'LineWidth', 2);
plot3(rrvvsEncRef(1), rrvvsEncRef(2), rrvvsEncRef(3), 'rs', 'LineWidth', 2);
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal
grid on
box on
savePlot(nameFig, path, type, res);

%% dv
figure(); hold on; grid on; box on
name = 'dvCompProf';
nameFig = strcat(tag,name);
yyaxis right
plot(tt,dvvRTN(1,1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', [0, 0.4470, 0.7410], 'linestyle', '-');
plot(tt,dvvRTN(2,1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], 'linestyle', '-');
plot(tt,dvvRTN(3,1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', [0.9290, 0.6940, 0.1250], 'linestyle', '-');
plot(tt,dv(1:end-1)*km/param.dt*mass*km, 'LineWidth', 2, 'color', 'k', 'linestyle', '-');
xlabel('Time to CA [rev]',  'color', 'k')
ylabel('Thrust [mN]',  'color', 'k')
set(gca,'xdir','reverse', 'ycolor', 'k')
legend('T', 'N', 'B', '\Delta v','location', 'southeast')
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
legend('T', 'N', 'B', '\Delta v','location', 'northwest')
ylim(1.1*[ylimmin ylimmax]./mass*param.dt)
savePlot(nameFig, path, type, res);

%% dv
figure(); hold on; grid on; box on
name = 'dvProf';
nameFig = strcat(tag,name);
semilogy(tt,dv(1:end-1)*km,'LineWidth', 2)
xlabel('Time to CA [rev]')
ylabel('dv [m/s]')
set(gca,'xdir','reverse')
savePlot(nameFig, path, type, res);