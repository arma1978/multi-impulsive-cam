function [out,xxOptLin,param] = linConvexSolveEllipse(param)
% linConvexSolveEllipse  Sweep a discretized B-plane ellipse and solve SOCP.
%
% For each sampled ellipse point, solve one convex problem and store
% delta-v and resulting B-plane state.
%
% Input:
%   param    - optimisation struct with ellipse samples and solver settings
% Outputs:
%   out      - struct array with per-sample optimisation results
%   xxOptLin - last optimisation decision vector (delta-v components)
%   param    - returned unchanged for interface consistency

for ind = 1:param.disc
    out(ind).xx2Bplane = param.xx2Target(:,ind);
    prob = linConvexProblemEllipse(param, out(ind).xx2Bplane);
    out(ind).xxOnEllipse = out(ind).xx2Bplane;
    tic
    [~,res]=mosekopt('minimize echo(0)',prob);
    toc
    assert(isfield(res, 'sol') && isfield(res.sol, 'itr'), ...
        'MOSEK returned no interior-point solution in linConvexSolveEllipse (rcode=%d: %s)', ...
        res.rcode, res.rcodestr);
    dum = res.sol.itr.xx;
    out(ind).dvv = reshape(dum(1:param.ndv*3), 3, param.ndv);
    out(ind).dv = dum(param.nOptVarDvComp+1:param.nOptVarDvComp+param.nOptVarDvMag);
    out(ind).dvTot = sum( out(ind).dv);
    out(ind).xx2Opt = dum(param.nOptVarDvComp+param.nOptVarDvMag+1:end);
    out(ind).tol = norm(out(ind).xx2Opt-out(ind).xxOnEllipse);
end

xxOptLin = dum(1:param.ndv*3);
