function [out, xxOptLin, DAx, param] = refConvex(dvv0, tcaLin, param)
% refConvex  One SCVX outer iteration: re-propagate around current solution,
%            rebuild DA maps, re-solve the linearised convex sub-problem.
%
% Inputs:
%   dvv0    - current best delta-v vector  [nOptVarDvComp x 1]  [km/s]
%   tcaLin  - current best TCA estimate   [s]
%   param   - parameter struct (updated in place)
%
% Outputs:
%   out      - struct array with per-inner-iteration data (B-plane points,
%              dv components, magnitudes, compute times)
%   xxOptLin - refined delta-v decision vector  [nOptVarDvComp x 1]  [km/s]
%   DAx      - updated DA flow map cell array (loaded from mapsRefine.dat)
%   param    - updated parameter struct

% Inner-loop convergence settings (B-plane point iteration)
tolRef    = 1e-3;   % convergence threshold on B-plane point movement [km]
maxInner  = 10;     % maximum inner iterations

%% Write propagation input file
% The binary reads tcaLin, the reference state, and the current delta-v
% sequence to re-propagate and rebuild DA maps centred on the new reference.
N     = min(floor(param.t2TCA / param.dt), param.maxImpNum);
nmaps = (N+2)*6 + 2 + 4 + 4;
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

fid = fopen(fullfile(runtimeDir,'xxs.dat'), 'w');
fprintf(fid, '%40.30f\n', tcaLin);
fprintf(fid, '%40.30f\n', param.xx0s(1));
fprintf(fid, '%40.30f\n', param.xx0s(2));
fprintf(fid, '%40.30f\n', param.xx0s(3));
fprintf(fid, '%40.30f\n', param.xx0s(4));
fprintf(fid, '%40.30f\n', param.xx0s(5));
fprintf(fid, '%40.30f\n', param.xx0s(6));
for kk = 1:param.nOptVarDvComp
    fprintf(fid, '%40.30f\n', dvv0(kk));
end
fclose(fid);

%% Run the propagator binary to produce mapsRefine.dat
if strcmp(param.prop, 'Kep')
    [propStatus, propOut] = system(fullfile(projectRoot,'bin','statePropMultiMapsKepRefine'));
    assert(propStatus == 0, 'Keplerian refinement propagator failed:\n%s', propOut);
else
    [propStatus, propOut] = system(fullfile(projectRoot,'bin','statePropMultiMapsFullPRefine'));
    assert(propStatus == 0, 'J2 refinement propagator failed:\n%s', propOut);
end

%% Load updated DA flow maps and rebuild linear quantities
param.fname = fullfile(runtimeDir,'mapsRefine.dat');
DAx = loadPoly(param.fname, param.nvar, nmaps, 0);

param = refTraj(DAx, param);       % reference trajectory at new expansion point
param = linearMaps(DAx, param);    % updated linear STMs
param = ellipseTargets(param);     % updated B-plane avoidance ellipse

%% Inner convex loop
%  Iteratively move the tangent point on the avoidance ellipse and re-solve
%  the SOCP until the B-plane point converges (SCVX inner loop).
out(1).xx2Bplane = pointOnEllipse(param.rrb2Nom, param);  % initial guess: closest point to nominal
ii = 1;  dist = Inf;
while (dist > tolRef) && (ii < maxInner)
    tic
    [prob, xxOnEllipse] = linConvexProblemRefine(param, out(ii).xx2Bplane, dvv0);
    out(ii).compTime   = toc;
    out(ii).xxOnEllipse = xxOnEllipse;

    % Solve the SOCP with MOSEK
    [~, res] = mosekopt('minimize echo(0)', prob);
    assert(isfield(res, 'sol') && isfield(res.sol, 'itr'), ...
        'MOSEK returned no interior-point solution in refConvex (rcode=%d: %s)', ...
        res.rcode, res.rcodestr);
    dum = res.sol.itr.xx;

    out(ii).dvv     = reshape(dum(1:param.ndv*3), 3, param.ndv);                          % [km/s] 3xN matrix
    out(ii).dv      = dum(param.nOptVarDvComp+1 : param.nOptVarDvComp+param.nOptVarDvMag); % magnitude aux vars
    out(ii).dvTot   = sum(out(ii).dv);                                                     % total dv [km/s]
    out(ii).xx2Opt  = dum(param.nOptVarDvComp+param.nOptVarDvMag+1 : end);                 % optimised B-plane point

    out(ii+1).xx2Bplane = out(ii).xx2Opt;
    dist = norm(out(ii).xx2Opt - xxOnEllipse);  % B-plane point movement this iteration
    ii   = ii + 1;
end
out(ii) = [];  % remove the uninitialised trailing entry

%% Extract final delta-v vector
xxOptLin = dum(1:param.ndv*3);
