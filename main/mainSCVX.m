% =========================================================================
% mainSCVX.m
% Sequential Convex (SCVX) Collision Avoidance Maneuver (CAM) optimisation
%
% Computes the minimum-delta-v impulsive manoeuvre sequence that drives the
% collision probability (or Mahalanobis distance) below a prescribed target,
% subject to a per-impulse thrust constraint.
%
% Workflow:
%   1.  Load conjunction catalogue and select / filter a test case
%   2.  Propagate the reference trajectory via an external C++ binary
%       and load the resulting differential-algebra (DA) flow maps
%   3.  Build and solve the linearised convex problem with MOSEK (1st round)
%   4.  Iteratively re-linearise and re-solve until convergence (SCVX loop)
%   5.  Post-process, print summary, generate figures, save workspace
%
% Units:  km, s  throughout (unless explicitly labelled otherwise)
% =========================================================================
clearvars; close all; clc; format longG;

% --- global figure formatting: LaTeX interpreter, Times font --------------
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultUicontrolFontName','Times', 'DefaultUicontrolFontSize', 14);
set(0,'DefaultUitableFontName','Times', 'DefaultUitableFontSize', 14);
set(0,'DefaultTextFontName','Times', 'DefaultTextFontSize', 14);
set(0,'DefaultUipanelFontName','Times', 'DefaultUipanelFontSize', 14);

% --- search paths ---------------------------------------------------------
projectRoot = fileparts(fileparts(mfilename('fullpath')));
cd(projectRoot);
addpath(fullfile(projectRoot,'src','routines'));
addpath(fullfile(projectRoot,'src','io'));
addpath(fullfile(projectRoot,'src','frames'));
addpath(fullfile(projectRoot,'utils'));
addpath(fullfile(getenv('HOME'),'Dropbox/Work/mosek/11.1/toolbox/r2022b')); % MOSEK 11 MATLAB toolbox
runtimeDir = fullfile(projectRoot,'runtime');
outputDir = fullfile(projectRoot,'output');
if ~exist(runtimeDir,'dir'), mkdir(runtimeDir); end

%% -------------------------------------------------------------------------
%  1. LOAD DATA
%     all units are km, s
% --------------------------------------------------------------------------
load(fullfile(projectRoot,'input','data.mat'))          % loads: data — conjunction catalogue (table)
% load dataRev.mat;    % alternative: larger simulation catalogue

%% -------------------------------------------------------------------------
%  2. PHYSICAL CONSTANTS AND SCALING
% --------------------------------------------------------------------------
km = 1000;                           % km-to-m factor (used only for display/output)
param.mu  = 398600.4418;             % Earth gravitational parameter  [km^3 s^-2]
param.rE  = 6378.137;                % Earth equatorial radius         [km]
param.J2  = 1.08262668e-3;           % J2 zonal harmonic coefficient   [-]
param.Lsc = param.rE;                % length scale for non-dimensionalisation [km]
param.Vsc = sqrt(param.mu/param.Lsc); % velocity scale                 [km/s]
param.Tsc = param.Lsc/param.Vsc;     % time scale                      [s]
param.disc = 1111;                   % number of discretisation points for DA map evaluation
mass = 300;                          % spacecraft wet mass             [kg]
Tmax = 10e-3;  % [N]  maximum continuous thrust (10 mN). Commented values show other test cases.

%% -------------------------------------------------------------------------
%  3. MISSION / OPTIMISATION PARAMETERS  (main user inputs)
% --------------------------------------------------------------------------
param.prop     = 'J2';  % dynamics model: 'J2' (oblate Earth) or 'Kep' (two-body)
param.nOrbits  = 8;     % manoeuvre horizon length [orbital periods of the chaser]
param.ind      = 1;     % row index of the selected conjunction in the catalogue
                        % (other tested cases: 53, 685, 889, 958, 1026, 1112)
param.order    = 2;     % order of the DA Taylor expansion
param.nvar     = 7;     % number of DA variables: 6 state + 1 time
param.dt       = 60;    % impulse spacing / propagation step size [s]
param.maxImpNum = 200;  % hard cap on number of impulses (avoids over-large problems)

% Maximum single-impulse delta-v derived from thrust level and step size
param.dvMax = Tmax / mass * param.dt / 1000;   % [km/s]  (N/kg * s / (m/km) )

% Risk thresholds used to select and constrain the optimisation
PcLim      = 1e-6;  % target Pc below which risk is acceptable (objective = 'risk')
maxPcLim   = 1e-4;  % maximum permissible Pc (hard constraint for other objectives)
minDistLim = 2;     % minimum acceptable miss distance  [km]

% Objective selector: 'miss_distance' | 'risk' | 'max_risk'
obj = 'max_risk';

% SCVX outer-loop convergence settings
tolSCVX      = 1e-6;  % sup-norm convergence threshold on the decision-vector correction
maxSCVXIter  = 10;    % hard cap on number of SCVX refinement iterations

%% -------------------------------------------------------------------------
%  4. OBJECTIVE FUNCTION AND FILTERING CONFIGURATION
%     Sets the risk/distance thresholds passed to selectTestCaseFiltered.
% --------------------------------------------------------------------------
switch obj
    case 'miss_distance'  % minimise miss distance; filter conjunctions below maxPcLim
        param.minDist = minDistLim;  % minimum target miss distance [km]
        order         = 'ascend';    % sort catalogue ascending by miss distance
        param.PcMax   = maxPcLim;
        param.Pc      = [];
    case 'risk'           % minimise collision probability; filter conjunctions above PcLim
        order         = 'descend';   % sort catalogue descending (highest risk first)
        param.Pc      = PcLim;       % target Pc
        param.minDist = [];
        param.PcMax   = [];
    case 'max_risk'       % minimise maximum Pc; filter conjunctions above maxPcLim
        order         = 'descend';   % sort catalogue descending (highest risk first)
        param.PcMax   = maxPcLim;
        param.minDist = [];
        param.Pc      = [];
    otherwise
        error('Unknown objective ''%s''. Choose: miss_distance | risk | max_risk', obj);
end

%% -------------------------------------------------------------------------
%  5. OUTPUT FILE PATHS
% --------------------------------------------------------------------------
param.outfile = fullfile(outputDir,'ResultsSCVX','sol');                     % prefix for saved result .mat files
param.fname   = fullfile(runtimeDir,'maps.dat');                             % binary file written by the propagator and read by loadPoly
param.projectRoot = projectRoot;
param.runtimeDir = runtimeDir;
if ~exist(fileparts(param.outfile),'dir'), mkdir(fileparts(param.outfile)); end

%% -------------------------------------------------------------------------
%  6. TEST CASE SELECTION
%     selectTestCaseFiltered filters the catalogue by the active risk/distance
%     thresholds, selects the conjunction at index param.ind, and populates
%     param with: initial orbital state (xx0), covariance matrices, and
%     time-to-closest-approach (t2TCA).
% --------------------------------------------------------------------------
param.data = data;
% param.data = sortrows(data,'d^* [km]', 'ascend');   % optional: pre-sort before filtering
param = selectTestCaseFiltered(param);

% Alternative: large-simulation data path (uses selectTestCase instead)
% param.data = dataRev;
% param.data = sortrows(param.data,'miss_distance', 'ascend');
% [param, covsRTN, covdRTN] = selectTestCase(param);

%% -------------------------------------------------------------------------
%  7. REFERENCE TRAJECTORY PROPAGATION  (external C++ binary)
%     The binary reads the initial conditions written by selectTestCaseFiltered,
%     propagates the spacecraft state and builds DA flow maps, then writes
%     results to disk (xs0*.dat and maps.dat).
%     param.xx0s: reference state at the start of the manoeuvre window [km; km/s]
% --------------------------------------------------------------------------
if strcmp(param.prop,'Kep')
    % Two-body Keplerian propagator + DA map builder
    [propStatus, propOut] = system(fullfile(projectRoot,'bin','statePropMultiMapsKep'));
    assert(propStatus == 0, 'Keplerian propagator failed:\n%s', propOut);
    dum = load(fullfile(runtimeDir,'xs0Kep.dat'));
else
    % J2-perturbed propagator + DA map builder
    [propStatus, propOut] = system(fullfile(projectRoot,'bin','statePropMultiMapsFullP'));
    assert(propStatus == 0, 'J2 propagator failed:\n%s', propOut);
    dum = load(fullfile(runtimeDir,'xs0J2.dat'));
end
param.xx0s = dum(1:6);   % [r (km); v (km/s)] at the beginning of the manoeuvre window
%% -------------------------------------------------------------------------
%  8. LOAD DA FLOW MAPS
%     N     = number of impulse opportunities (time steps until TCA, capped at maxImpNum)
%     nmaps = total number of map columns in the binary file:
%               N*6  state maps (6 state components per step)
%               2*6  maps at the encounter epoch (chaser + debris)
%               2    TCA scalar  +  Mahalanobis distance scalar
%               4    Mahalanobis matrix elements
%               4    covariance-related terms
%     DAx   = cell array of DA polynomial maps loaded from maps.dat
% --------------------------------------------------------------------------
N     = min(floor(param.t2TCA / param.dt), param.maxImpNum);
nmaps = (N+2)*6 + 2 + 4 + 4;
DAx   = loadPoly(param.fname, param.nvar, nmaps, 0);

%% -------------------------------------------------------------------------
%  9. OPTIMISATION VARIABLE LAYOUT
%     The MOSEK decision vector is structured as:
%       x = [ dv_xyz (N*3) | dv_mag (N) | state_corrections (2) ]
%     dv_xyz:             3 Cartesian components per impulse node
%     dv_mag:             1 auxiliary magnitude variable per node (used in SOC constraints)
%     state_corrections:  2 scalars encoding the TCA and Mahalanobis correction
% --------------------------------------------------------------------------
nOptVarDvComp       = N * 3;   % number of delta-v component variables
nOptVarDvMag        = N;       % number of delta-v magnitude auxiliary variables
nOptVarState        = 2;       % number of state correction variables
nOptVar             = nOptVarDvComp + nOptVarDvMag + nOptVarState;  % total decision variables
param.nOptVar       = nOptVar;
param.nOptVarDvComp = nOptVarDvComp;
param.nOptVarDvMag  = nOptVarDvMag;
param.nOptVarState  = nOptVarState;
param.ndv           = nOptVarDvComp / 3;   % number of impulse nodes = N

%% -------------------------------------------------------------------------
%  10. PRE-COMPUTATION
%      refTraj:       evaluates DA maps at the reference point → reference
%                     state and velocity history along the manoeuvre window
%      linearMaps:    extracts the first-order (linear) state-transition matrices
%                     from the DA maps (Jacobians w.r.t. impulse positions)
%      ellipseTargets: constructs the collision-avoidance target ellipse in the
%                     B-plane from the covariance matrices and Pc threshold
%      validate:      computes risk metrics for a given delta-v sequence; called
%                     here with zeros to get baseline (uncontrolled) metrics
% --------------------------------------------------------------------------
param = refTraj(DAx, param);       % reference trajectory positions and velocities
param = linearMaps(DAx, param);    % linear STMs from each impulse node to TCA
param = ellipseTargets(param);     % B-plane avoidance ellipse geometry

% Baseline metrics (zero manoeuvre): establishes the pre-CAM risk level
[sqrMaha0, tca0, dv0, dvv0, xx0, xxNoMan0, xxsEnc0, xxdEnc0, ...
    rrb2d0, Pb_2D0, maxPc0, Pc0] = validate(zeros(nOptVarDvComp, 1), DAx, param);

%% -------------------------------------------------------------------------
%  11. FIRST CONVEX SOLVE  (linearised problem, solved by MOSEK)
%      linConvexSolve builds the SOCP with thrust-magnitude SOC constraints,
%      a linear Mahalanobis distance constraint, and a minimum total Pc
%      objective, then calls mosekopt.
%      outLin:   struct array (one entry per MOSEK iteration) with
%                intermediate solutions and convergence data
%      xxOptLin: optimal decision vector  [length nOptVar]
% --------------------------------------------------------------------------
tic
[outLin, xxOptLin, param] = linConvexSolve(param);
compTime = toc;
fprintf('1st-round solve time: %.4f s\n', compTime);

% Validate the 1st-round solution against the nonlinear dynamics
[sqrMahalanobisLin, tcaLin, dvLin, dvvLin, ...
    rrvvLin, rrvvNoManLin, rrvvsEncLin, rrvvdEncLin, ...
    rrb2dLin, Pb_2DLin, maxPcLin, PcLin] = validate(xxOptLin, DAx, param);

%% -------------------------------------------------------------------------
%  12. SCVX REFINEMENT LOOP
%      Each iteration re-propagates around the current solution (refConvex),
%      re-linearises the DA maps, and re-solves the convex sub-problem.
%      The increment xx(:,jj+1) - xx(:,jj) represents the correction added
%      in each SCVX step; validate receives this increment.
%      Convergence criterion: max element-wise change in the decision vector
%      tolM < 1e-6, or a hard cap of 10 outer iterations.
% --------------------------------------------------------------------------
tolM    = Inf;         % initialise above threshold to enter the loop unconditionally
xx(:,1) = xxOptLin;   % column 1: solution from 1st-round convex solve
tca(1)  = tcaLin;
pp(1)   = param;
jj = 0;
while (tolM > tolSCVX) && (jj < maxSCVXIter)
    jj = jj + 1;
    % Re-linearise around current solution and solve the updated convex problem
    [outRef, xx(:,jj+1), DAxRef, pp(jj+1)] = refConvex(xx(:,jj), tca(jj), pp(jj));
    compTime = compTime + sum([outRef.compTime]);
    % validate expects the step increment (correction), not the full state
    [~, tca(jj+1), ~, ~, ~, ~, ~, ~, rrb2dRefIter(:,jj)] = ...
        validate(xx(:,jj+1) - xx(:,jj), DAxRef, pp(jj+1));
    tolM = max(abs(xx(:,jj+1) - xx(:,jj)));   % sup-norm of the correction
end
paramRef = pp(end);   % param struct corresponding to the final linearisation point

%% -------------------------------------------------------------------------
%  13. POST-PROCESSING  (final refined solution)
% --------------------------------------------------------------------------
% Validate the final SCVX correction step against the nonlinear dynamics
[sqrMahalanobisRef, tcaRef, dvRef, dvvRef, ...
    rrvvRef, rrvvNoManRef, rrvvsEncRef, rrvvdEncRef, ...
    rrb2dRef, Pb_2DRef, maxPcRef, PcRef] = ...
    validate(xx(:,jj+1) - xx(:,jj), DAxRef, paramRef);

% Reconstruct full impulse matrix (3 x N+1); append a zero impulse at TCA
dvv = [reshape(xx(:,end), 3, param.ndv)  zeros(3,1)];   % [km/s] each column
dv  = sqrt(sum(dvv .* dvv, 1));                         % impulse magnitudes [km/s]

% Rotation matrices from inertial to Frenet (RTN) frame at each impulse node
% frenet() is preferred over RTN() for non-circular reference orbits
toRTN = frenet(rrvvRef);

% Project impulses into RTN frame and compute classical orbital elements
% Pre-allocate outputs for clarity (loop over all nodes including the zero TCA one)
dvvRTN = zeros(3, size(dvv,2));
coe    = zeros(6, size(dvv,2));
for kk = 1:size(dvv, 2)
    dvvRTN(:,kk) = toRTN(:,:,kk) * dvv(:,kk);                           % [km/s] in RTN
    coe(:,kk)    = pv2po(rrvvRef(1:3,kk), rrvvRef(4:6,kk), param.mu);   % [km, rad, ...]
end

% Build time-to-CA axis in orbital revolutions for the delta-v plots
T  = 2*pi * sqrt(coe(1, 1:end-1).^3 / param.mu);   % local orbital period at each node [s]
tt = param.nOrbits - cumsum(param.dt ./ T);          % time to CA  [revolutions]

%% -------------------------------------------------------------------------
%  14. OUTPUT SUMMARY
% --------------------------------------------------------------------------
fprintf('\n--- SCVX / MOSEK convergence summary ---\n');
fprintf('  Solve time (total):        %.4f s\n', compTime);
fprintf('  MOSEK iterations:          %i (1st round),  %i (refinement round %i)\n', ...
    size(outLin,2), size(outRef,2), jj+1);
fprintf('  Total delta-v:             %.4f mm/s  (1st round)\n', sum(dvLin)*km^2);
fprintf('                             %.4f mm/s  (round %i)\n', sum(dv)*km^2, jj+1);
fprintf('  Squared Mahalanobis dist:  %.4f  (1st round)\n', sqrMahalanobisLin);
fprintf('                             %.4f  (round %i)\n', sqrMahalanobisRef, jj+1);
fprintf('  Target Mahalanobis dist:   %.4f\n', param.sqrMahaMin);
fprintf('  Initial TCA:               %.4f s\n', param.t2TCA);
fprintf('  TCA after 1st round:       %.4f s\n', tcaLin);
fprintf('  TCA after round %i:         %.4f s\n', jj+1, tcaRef);

%% -------------------------------------------------------------------------
%  15. PLOTS
% --------------------------------------------------------------------------

% --- Fig 1: B-plane view
%     Shows the optimised solution positions projected onto the collision
%     avoidance ellipse in the B-plane (zeta, xi coordinates).
%     Black/blue curves: target avoidance ellipses (1st round vs. refined).
%     Scatter: impulse nodes, coloured by total delta-v [m/s].
%     Red circles: closest point on ellipse for each iteration.
%     Black/green dots: nominal and achieved conjunction B-plane point.
figure(); hold on; grid on; box on
plot(param.xx2Target(2,:),    param.xx2Target(1,:),    'k', 'LineWidth', 2)  % 1st-round ellipse
plot(paramRef.xx2Target(2,:), paramRef.xx2Target(1,:), 'b', 'LineWidth', 2)  % refined ellipse
dumT = [outLin.xx2Opt  outRef.xx2Opt];            % optimised positions in B-plane
dumv = [outLin.dvTot   outRef.dvTot] * 1000;      % total dv per iteration [m/s]
for ii = 1:length(dumv),  textstr{ii} = num2str(ii);  end
text(dumT(2,:), dumT(1,:), textstr, 'FontSize', 20, 'FontName', 'Times')
scatter(dumT(2,:), dumT(1,:), 200, dumv', 'filled', 'MarkerEdgeColor', 'k')
dumT = [outLin.xxOnEllipse  outRef.xxOnEllipse];   % nearest points on ellipse
plot(dumT(2,:), dumT(1,:), 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5)
text(dumT(2,:), dumT(1,:), textstr, 'FontSize', 20, 'FontName', 'Times')
plot(param.rrb2Nom(2), param.rrb2Nom(1), 'ko', 'MarkerFaceColor','k', 'MarkerEdgeColor','k')  % nominal (no-manoeuvre) point
plot(rrb2dRef(2),      rrb2dRef(1),      'ko', 'MarkerFaceColor','g', 'MarkerEdgeColor','g')  % achieved point
ylabel('$\xi$ [km]');  xlabel('$\zeta$ [km]')
colorbar;  title('$\Delta v$ [m/s]');  axis equal

% --- Fig 2: 3-D inertial trajectory with impulse quiver arrows
%     Red arrows: delta-v impulses at each node (quiver3, not to scale).
%     Blue curve:  uncontrolled (no-manoeuvre) trajectory.
%     Red square:  chaser position at closest approach (TCA).
figure(); hold on; grid on; box on
quiver3(rrvvRef(1,:), rrvvRef(2,:), rrvvRef(3,:), ...
        dvv(1,:),    dvv(2,:),    dvv(3,:), 'r', 'LineWidth', 2);
plot3(rrvvNoManRef(1,1:end-1), rrvvNoManRef(2,1:end-1), rrvvNoManRef(3,1:end-1), ...
      'b', 'LineWidth', 2);
plot3(rrvvsEncRef(1), rrvvsEncRef(2), rrvvsEncRef(3), 'rs', 'LineWidth', 2);
xlabel('x [km]');  ylabel('y [km]');  zlabel('z [km]')
axis equal;  grid on;  box on

% --- Fig 3: RTN delta-v components vs time to CA  [m/s, linear scale]
figure(); hold on; grid on; box on
plot(tt, dvvRTN(1,1:end-1)*km, 'LineWidth', 2);   % Tangential
plot(tt, dvvRTN(2,1:end-1)*km, 'LineWidth', 2);   % Normal
plot(tt, dvvRTN(3,1:end-1)*km, 'LineWidth', 2);   % Bi-normal
plot(tt, dv(1:end-1)*km,       'LineWidth', 2);   % Total magnitude
xlabel('Time to CA [rev]');  ylabel('$\Delta v$ [m/s]')
set(gca,'xdir','reverse');   legend('T', 'N', 'B', '\Delta v')

% --- Fig 4: delta-v magnitude vs time to CA  [m/s, log scale]
figure(); hold on; grid on; box on
semilogy(tt, dv(1:end-1)*km, 'LineWidth', 2)
xlabel('Time to CA [rev]');  ylabel('$\Delta v$ [m/s]')
set(gca,'xdir','reverse')

% --- Fig 5: Dual-axis — thrust [mN] (right axis) and delta-v [mm/s] (left axis)
%     Conversion:  thrust [mN] = (dv [km/s] / dt [s]) * mass [kg] * 1e6
%                  dv    [mm/s] = dv [km/s] * km^2  (km^2 = 1e6)
figure(); hold on; grid on; box on

% Shared colour scheme for T / N / B / total
cT = [0, 0.4470, 0.7410];    % blue  — Tangential
cN = [0.8500, 0.3250, 0.0980]; % orange — Normal
cB = [0.9290, 0.6940, 0.1250]; % yellow — Bi-normal

% Common scaling factor: dv [km/s] → thrust [mN]
thrustScale = km / param.dt * mass * km;  % (m/km) / s * kg * (m/km) = N * 1e6 = mN... check units

yyaxis right
plot(tt, dvvRTN(1,1:end-1)*thrustScale, 'LineWidth', 2, 'color', cT, 'linestyle', '-');
plot(tt, dvvRTN(2,1:end-1)*thrustScale, 'LineWidth', 2, 'color', cN, 'linestyle', '-');
plot(tt, dvvRTN(3,1:end-1)*thrustScale, 'LineWidth', 2, 'color', cB, 'linestyle', '-');
plot(tt, dv(1:end-1)*thrustScale,        'LineWidth', 2, 'color', 'k', 'linestyle', '-');
xlabel('Time to CA [rev]', 'color', 'k');  ylabel('Thrust [mN]', 'color', 'k')
set(gca, 'xdir', 'reverse', 'ycolor', 'k')
legend('T', 'N', 'B', '\Delta v', 'location', 'northwest')

% Compute axis limits from all RTN components together
allThrust = [dvvRTN(1,1:end-1)  dvvRTN(2,1:end-1)  dvvRTN(3,1:end-1)  dv(1:end-1)] * thrustScale;
ylimmax = max(allThrust);  ylimmin = min(allThrust);
ylim(1.1 * [ylimmin  ylimmax])

yyaxis left
plot(tt, dvvRTN(1,1:end-1)*km^2, 'LineWidth', 2, 'color', cT, 'linestyle', '-');
plot(tt, dvvRTN(2,1:end-1)*km^2, 'LineWidth', 2, 'color', cN, 'linestyle', '-');
plot(tt, dvvRTN(3,1:end-1)*km^2, 'LineWidth', 2, 'color', cB, 'linestyle', '-');
plot(tt, dv(1:end-1)*km^2,        'LineWidth', 2, 'color', 'k', 'linestyle', '-');
xlabel('Time to CA [rev]', 'color', 'k');  ylabel('$\Delta v$ [mm/s]', 'color', 'k')
set(gca, 'xdir', 'reverse', 'ycolor', 'k')
legend('T', 'N', 'B', '\Delta v', 'location', 'northeast')
ylim(1.1 * [ylimmin  ylimmax] ./ mass * param.dt)


%% -------------------------------------------------------------------------
%  16. SAVE WORKSPACE
%     Large DA arrays (data, DAx, DAxRef) are cleared before saving to reduce
%     file size; they can be regenerated from param settings if needed.
% --------------------------------------------------------------------------
clear data DAx DAxRef
fnamemat = strcat([param.outfile, param.prop, '_ID_', num2str(param.ind), ...
    '_obj_', obj, '_time2man_', num2str(param.nOrbits), '_maxnImp_', num2str(param.maxImpNum)]);
save(fnamemat);