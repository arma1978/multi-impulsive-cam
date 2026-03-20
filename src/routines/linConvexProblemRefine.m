function [prob, xxOnEllipse] = linConvexProblemRefine(param, xx0Guess, dvv0)
% linConvexProblemRefine  Build the linearised SOCP for SCVX refinement iterations.
%
%   Same structure as linConvexProblem, but the constraints are expressed
%   as corrections relative to the current linearisation point dvv0.
%   The decision variable is therefore delta_dv = dv - dvv0.
%
%   Inequality (corrected for dvv0 offset):
%     -n' * drrb2Nom * ddv <= n' * (rrb2Nom - xxOnEllipse - drrb2Nom*dvv0)
%
%   Equality (corrected for dvv0 offset):
%     drrb2Nom * ddv - bplane_pt = -rrb2Nom + drrb2Nom*dvv0
%
% Inputs:
%   param     - struct with updated linear maps at new expansion point
%   xx0Guess  - current B-plane point guess  [2 x 1]  [km]
%   dvv0      - current reference delta-v vector  [nOptVarDvComp x 1]  [km/s]
% Outputs:
%   prob        - MOSEK problem struct
%   xxOnEllipse - nearest point on avoidance ellipse to xx0Guess  [2 x 1]  [km]

dvmax = param.dvMax;
[xxOnEllipse, xxOnEllipseNormal] = pointOnEllipse(xx0Guess, param);

%% Inequality: tangent half-plane constraint (corrected for linearisation offset)
A = -xxOnEllipseNormal' * param.drrb2Nom;
A = [A  zeros(1, param.nOptVarDvMag)  zeros(1, param.nOptVarState)];
b =  xxOnEllipseNormal' * (param.rrb2Nom - xxOnEllipse - param.drrb2Nom*dvv0);

%% Equality: B-plane state variable (corrected for linearisation offset)
Aeq = [param.drrb2Nom  zeros(2, param.nOptVarDvMag)  -eye(2, param.nOptVarState)];
beq = -param.rrb2Nom + param.drrb2Nom*dvv0;

%% Assemble MOSEK problem struct
prob.c   = [zeros(1, param.nOptVarDvComp)  ones(1, param.nOptVarDvMag)  zeros(1, param.nOptVarState)];
prob.a   = sparse([A; Aeq]);
prob.blc = [-inf; beq];
prob.buc = [b;    beq];
prob.blx = -[dvmax * ones(1, param.nOptVarDvComp + param.nOptVarDvMag)  200  200];
prob.bux = -prob.blx;

%% Second-order cone constraints  (identical layout to linConvexProblem)
prob.cones.type   = zeros(1, param.nOptVarDvMag);
indices  = 1 : param.nOptVarDvComp + param.nOptVarDvMag;
dvindec  = reshape(indices(1:param.nOptVarDvComp), 3, param.nOptVarDvComp/3)';
dvmag    = indices(param.nOptVarDvComp+1:end)';
prob.cones.sub    = reshape([dvmag  dvindec]', param.nOptVarDvComp + param.nOptVarDvMag, 1)';
prob.cones.subptr = 1 : 4 : param.nOptVarDvComp + param.nOptVarDvMag;

end