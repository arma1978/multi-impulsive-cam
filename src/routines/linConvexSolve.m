function [out, xxOptLin, param] = linConvexSolve(param)
% linConvexSolve  First-round linearised convex solve.
%
%   Iterates the B-plane tangent-point procedure: at each step the nearest
%   point on the avoidance ellipse is recomputed and the SOCP is re-solved
%   until the optimised B-plane point converges.  Uses the DA maps already
%   stored in param (computed once in mainSCVX before this call).
%
% Input:
%   param    - parameter struct including linear maps, ellipse geometry,
%              optimisation variable sizes, and dvMax
% Outputs:
%   out      - struct array (one entry per iteration) with B-plane points,
%              dv components, magnitudes
%   xxOptLin - optimal delta-v decision vector  [nOptVarDvComp x 1]  [km/s]
%   param    - parameter struct (unchanged by this function)

% Inner-loop convergence settings
tolInner  = 1e-3;   % B-plane point convergence threshold [km]
maxInner  = 10;     % maximum inner iterations

%% B-plane tangent-point iteration
%  Initialise at the point on the avoidance ellipse closest to the nominal
%  (uncontrolled) conjunction point.
out(1).xx2Bplane = pointOnEllipse(param.rrb2Nom, param);

ii = 1;  tolm = Inf;
while (tolm > tolInner) && (ii < maxInner)
    % Build the linearised SOCP around the current B-plane tangent point
    [prob, xxOnEllipse]  = linConvexProblem(param, out(ii).xx2Bplane);
    out(ii).xxOnEllipse  = xxOnEllipse;

    % Solve with MOSEK
    [~, res] = mosekopt('minimize echo(0)', prob);
    assert(isfield(res, 'sol') && isfield(res.sol, 'itr'), ...
        'MOSEK returned no interior-point solution in linConvexSolve (rcode=%d: %s)', ...
        res.rcode, res.rcodestr);
    dum = res.sol.itr.xx;

    out(ii).dvv    = reshape(dum(1:param.ndv*3), 3, param.ndv);                           % [km/s] 3xN
    out(ii).dv     = dum(param.nOptVarDvComp+1 : param.nOptVarDvComp+param.nOptVarDvMag); % magnitude aux vars
    out(ii).dvTot  = sum(out(ii).dv);                                                     % total dv [km/s]
    out(ii).xx2Opt = dum(param.nOptVarDvComp+param.nOptVarDvMag+1 : end);                 % B-plane point

    out(ii+1).xx2Bplane = out(ii).xx2Opt;
    tolm = norm(out(ii).xx2Opt - xxOnEllipse);  % B-plane point movement this iteration
    ii   = ii + 1;
end
out(ii) = [];  % remove the uninitialised trailing entry

%% Extract final delta-v vector
xxOptLin = dum(1:param.ndv*3);
