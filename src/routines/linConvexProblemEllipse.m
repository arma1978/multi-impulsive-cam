function prob = linConvexProblemEllipse(param, xxOnellipse)
% linConvexProblemEllipse  Build SOCP with B-plane point fixed on ellipse grid.
%
% Inputs: %   param        - optimisation and linear-map data
%   xxOnellipse  - fixed B-plane target point [2x1]
% Outputs: %   prob         - MOSEK problem structure

dvmax = param.dvMax;

% Equality constraints (2): map delta-v to B-plane state variables
Aeq = [param.drrb2Nom zeros(2,param.nOptVarDvMag) -eye(2, param.nOptVarState)];
beq = -param.rrb2Nom;

% See: https://docs.mosek.com/latest/toolbox/tutorial-cqo-shared.html

% Specify the non-conic part of the problem.
prob.c   = [zeros(1, param.nOptVarDvComp) ones(1,param.nOptVarDvMag) zeros(1, param.nOptVarState)]; % function to be minimised

prob.a   = sparse(Aeq);
prob.blc = beq;
prob.buc = beq;
prob.blx = [-dvmax*ones(1,param.nOptVarDvComp+param.nOptVarDvMag) xxOnellipse(1,1) xxOnellipse(2,1)];
prob.bux = [dvmax*ones(1,param.nOptVarDvComp+param.nOptVarDvMag)  xxOnellipse(1,1) xxOnellipse(2,1)];
% Specify the cones.    
prob.cones.type   = zeros(1,param.nOptVarDvMag);
indices = 1:param.nOptVarDvComp + param.nOptVarDvMag;
dvindex = reshape(indices(1:param.nOptVarDvComp), 3, param.nOptVarDvComp/3)';
dvmag   = indices(param.nOptVarDvComp+1:end)';
prob.cones.sub    = reshape([dvmag dvindex]',param.nOptVarDvComp + param.nOptVarDvMag,1)';
prob.cones.subptr = 1:4:param.nOptVarDvComp + param.nOptVarDvMag; % starting index of each cone

end

