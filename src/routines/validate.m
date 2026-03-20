function [sqrMahalanobis, tca, dv, dvv, ...
    rrvv, rrvvNoMan, rrvvsEnc, rrvvdEnc, ...
    rrb2d, Pb_2D, maxPc, Pc] = validate(xx, DAx, param)
% validate  Evaluate risk metrics for a given delta-v sequence using DA maps.
%
%   Applies the delta-v impulses to the DA flow maps, evaluates the
%   resulting encounter state, B-plane position, covariance, and computes
%   both the squared Mahalanobis distance and collision probability.
%
% Inputs:
%   xx    - delta-v component vector  [nOptVarDvComp x 1]  [km/s]
%           (pass zeros to evaluate the uncontrolled baseline)
%   DAx   - DA flow map cell array (from loadPoly)
%   param - parameter struct
%
% Outputs:
%   sqrMahalanobis - squared Mahalanobis distance at TCA  [-]
%   tca            - time of closest approach             [s]
%   dv             - impulse magnitudes  [1 x N+1]         [km/s]
%   dvv            - impulse vectors     [3 x N+1]         [km/s]
%   rrvv           - chaser state history (controlled)     [6 x N+1]  [km; km/s]
%   rrvvNoMan      - chaser state history (unmanoeuvred)   [6 x N+1]
%   rrvvsEnc       - chaser state at encounter             [6 x 1]
%   rrvvdEnc       - debris  state at encounter             [6 x 1]
%   rrb2d          - B-plane displacement vector           [2 x 1]    [km]
%   Pb_2D          - 2D B-plane covariance matrix          [2 x 2]    [km^2]
%   maxPc          - maximum Pc (upper bound)              [-]
%   Pc             - collision probability                 [-]

%% Unpack decision vector
ndv = param.ndv;
dvv = reshape(xx(1:ndv*3), 3, ndv);           % [km/s]  3 x N impulse matrix
DV  = [zeros(size(dvv)); dvv; zeros(1,ndv)];  % 7 x N  (pad with zero time variation)
DAMaps = reshape(DAx(1:ndv*6), 6, ndv);       % 6 x N  DA polynomial maps

% Pre-allocate state arrays
DX        = zeros(7, ndv);
dxf       = zeros(6, 1);
dxfNoMan  = dxf;
rrvvsEnc  = dxf;
rrvvdEnc  = dxf;
rrvv      = zeros(6, ndv);
rrvvNoMan = rrvv;

%% Propagate state through impulse sequence using DA maps
%  DX accumulates the state deviation due to the delta-v sequence.
%  At each node the map is evaluated at dx0 = accumulated deviation + impulse.
for ii = 1:ndv
    dx0      = DX(:,ii)' + DV(:,ii)';   % perturbed state including current impulse
    dx0NoMan = zeros(1,7);               % unperturbed reference (no manoeuvre)
    for jj = 1:6
        dxf(jj)      = evalPoly(DAMaps(jj,ii).C, DAMaps(jj,ii).E, dx0);
        dxfNoMan(jj) = evalPoly(DAMaps(jj,ii).C, DAMaps(jj,ii).E, dx0NoMan);
    end
    DX(1:6, ii+1) = dxf - dxfNoMan;  % incremental deviation from no-manoeuvre trajectory
    rrvv(:,ii)      = dxf;
    rrvvNoMan(:,ii) = dxfNoMan;
end

%% Evaluate encounter (TCA) state from the propagated final deviation
dx0         = DX(:,end)';
DAEncounter = reshape(DAx(ndv*6+1 : ndv*6+12), 6, 2);  % 6 x 2: [chaser, debris] at encounter
for jj = 1:6
    rrvvsEnc(jj,1) = evalPoly(DAEncounter(jj,1).C, DAMaps(jj,1).E, dx0);  % chaser at TCA
    rrvvdEnc(jj,1) = evalPoly(DAEncounter(jj,2).C, DAMaps(jj,2).E, dx0);  % debris at TCA
end

%% Compute actual TCA (DA map captures nonlinear TCA variation)
DAtca = DAx(ndv*6+12+1);
tca   = evalPoly(DAtca.C, DAtca.E, dx0);

%% Compute B-plane covariance at actual encounter
%  Indices: +7 = TCA, +8..+11 = Pb_2D elements (11=off-diag, symmetrised)
DAPb_2D    = DAx(ndv*6+12+7 : ndv*6+12+10);
Pb_2D(1,1) = evalPoly(DAPb_2D(1).C, DAPb_2D(1).E, dx0);
Pb_2D(1,2) = evalPoly(DAPb_2D(2).C, DAPb_2D(2).E, dx0);
Pb_2D(2,1) = Pb_2D(1,2);
Pb_2D(2,2) = evalPoly(DAPb_2D(4).C, DAPb_2D(4).E, dx0);

%% Compute squared Mahalanobis distance in the B-plane
%  rrb2d      = relative position in B-plane  [km]
%  Nb2drrb2d  = Nb_2D * rrb2d  (covariance-weighted direction)
DArrb2d     = DAx(ndv*6+12+3 : ndv*6+12+4);
DANb2drrb2d = DAx(ndv*6+12+5 : ndv*6+12+6);

rrb2d(1,1)     = evalPoly(DArrb2d(1).C,     DArrb2d(1).E,     DX(:,end)');
rrb2d(2,1)     = evalPoly(DArrb2d(2).C,     DArrb2d(2).E,     DX(:,end)');
Nb2drrb2d(1,1) = evalPoly(DANb2drrb2d(1).C, DANb2drrb2d(2).E, DX(:,end)');
Nb2drrb2d(2,1) = evalPoly(DANb2drrb2d(2).C, DANb2drrb2d(2).E, DX(:,end)');

sqrMahalanobis = rrb2d' * Nb2drrb2d;

%% Compute collision probability metrics
%  maxPc: maximum Pc upper bound (Chan formula)  [-]
%  Pc:    Gaussian collision probability          [-]
maxPc = param.bodySize^2 / exp(1) / sqrMahalanobis / det(Pb_2D)^0.5;
Pc    = param.bodySize^2 / 2      / det(Pb_2D)^0.5  * exp(-sqrMahalanobis/2);

%% Append zero impulse at TCA and prepend initial state
dvv       = [dvv  zeros(3,1)];          % [3 x N+1]
rrvv      = [param.xx0s  rrvv];         % [6 x N+1]
rrvvNoMan = [param.xx0s  rrvvNoMan];    % [6 x N+1]
dv        = sqrt(sum(dvv .* dvv, 1));   % impulse magnitudes  [1 x N+1]
end