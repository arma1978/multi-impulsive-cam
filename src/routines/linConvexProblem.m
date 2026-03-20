function [prob, xxOnEllipse] = linConvexProblem(param, xx0Guess)
% linConvexProblem  Build the linearised SOCP for the first-round convex solve.
%
%   Decision vector layout:  x = [ dv_xyz (nDvComp) | dv_mag (nDvMag) | bplane_pt (2) ]
%
%   Objective:  minimise sum of dv_mag  (total delta-v)
%
%   Constraints:
%     Inequality: B-plane point must lie in the half-plane tangent to the
%                 avoidance ellipse at xxOnEllipse (linearised ellipse constraint)
%     Equality:   B-plane point = drrb2Nom * dv_xyz + rrb2Nom
%                 (linear map from delta-v to B-plane displacement)
%     Box:        |dv_xyz| <= dvMax (per-component),  |bplane_pt| <= 200 km
%     SOC:        ||dv_xyz(3k-2:3k)|| <= dv_mag(k)  for each impulse k
%                 (second-order cone linking components to magnitude)
%
% See: https://docs.mosek.com/latest/toolbox/tutorial-cqo-shared.html
%
% Inputs:
%   param     - struct with drrb2Nom, rrb2Nom, dvMax, nOptVar* fields
%   xx0Guess  - initial guess for the B-plane point  [2 x 1]  [km]
% Outputs:
%   prob        - MOSEK problem struct ready for mosekopt
%   xxOnEllipse - nearest point on avoidance ellipse to xx0Guess  [2 x 1]  [km]

dvmax = param.dvMax;
[xxOnEllipse, xxOnEllipseNormal] = pointOnEllipse(xx0Guess, param);

%% Inequality: B-plane point in tangent half-plane at xxOnEllipse
%  n' * (drrb2Nom*dv + rrb2Nom) >= n' * xxOnEllipse
%  (written as -n'*drrb2Nom*dv <= n'*(rrb2Nom - xxOnEllipse))
A = -xxOnEllipseNormal' * param.drrb2Nom;
A = [A  zeros(1, param.nOptVarDvMag)  zeros(1, param.nOptVarState)];
b =  xxOnEllipseNormal' * (param.rrb2Nom - xxOnEllipse);

%% Equality: B-plane state variable equals linear map of delta-v
%  drrb2Nom * dv_xyz - bplane_pt = -rrb2Nom
Aeq = [param.drrb2Nom  zeros(2, param.nOptVarDvMag)  -eye(2, param.nOptVarState)];
beq = -param.rrb2Nom;

%% Assemble MOSEK problem struct
prob.c   = [zeros(1, param.nOptVarDvComp)  ones(1, param.nOptVarDvMag)  zeros(1, param.nOptVarState)];
prob.a   = sparse([A; Aeq]);
prob.blc = [-inf; beq];
prob.buc = [b;    beq];
prob.blx = -[dvmax * ones(1, param.nOptVarDvComp + param.nOptVarDvMag)  200  200];
prob.bux = -prob.blx;

%% Second-order cone constraints:  ||dv_xyz_k|| <= dv_mag_k  for each impulse k
%  MOSEK cone format: [mag; x; y; z] per cone block
prob.cones.type   = zeros(1, param.nOptVarDvMag);   % 0 = quadratic cone
indices  = 1 : param.nOptVarDvComp + param.nOptVarDvMag;
dvindec  = reshape(indices(1:param.nOptVarDvComp), 3, param.nOptVarDvComp/3)';  % N x 3
dvmag    = indices(param.nOptVarDvComp+1:end)';                                 % N x 1
prob.cones.sub    = reshape([dvmag  dvindec]', param.nOptVarDvComp + param.nOptVarDvMag, 1)';
prob.cones.subptr = 1 : 4 : param.nOptVarDvComp + param.nOptVarDvMag;  % start index of each cone block

end

