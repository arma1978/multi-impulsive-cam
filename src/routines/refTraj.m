function param = refTraj(DAx, param)
% refTraj  Evaluate the unperturbed (no-manoeuvre) reference trajectory
%          at the zero-deviation point of the DA maps.
%
%   Evaluates each DA flow map at dx0 = 0 to obtain the reference state
%   at every impulse node and the encounter states at TCA.
%   Results are stored in param for use by linearMaps, ellipseTargets,
%   and validate.
%
% Inputs:
%   DAx   - DA flow map cell array (from loadPoly)
%   param - struct with field ndv
% Output:
%   param - updated with rrvvNoMan, rrvvsEncNoMan, rrvvdEncNoMan

ndv     = param.ndv;
DAMaps  = reshape(DAx(1:ndv*6), 6, ndv);  % 6 x N DA maps along the manoeuvre window
dx0NoMan = zeros(1, 7);                    % zero-deviation expansion point

%% Reference trajectory at each impulse node  [km; km/s]
rrvvNoMan = zeros(6, ndv);
for ii = 1:ndv
    for jj = 1:6
        rrvvNoMan(jj,ii) = evalPoly(DAMaps(jj,ii).C, DAMaps(jj,ii).E, dx0NoMan);
    end
end
param.rrvvNoMan = rrvvNoMan;

%% Reference encounter states (chaser and debris) at TCA
DAEncounter = reshape(DAx(ndv*6+1 : ndv*6+12), 6, 2);  % columns: [chaser | debris]
dx0 = zeros(1, 7);
for jj = 1:6
    rrvvsEnc(jj,1) = evalPoly(DAEncounter(jj,1).C, DAEncounter(jj,1).E, dx0);  % chaser
    rrvvdEnc(jj,1) = evalPoly(DAEncounter(jj,2).C, DAEncounter(jj,2).E, dx0);  % debris
end
param.rrvvsEncNoMan = rrvvsEnc;
param.rrvvdEncNoMan = rrvvdEnc;