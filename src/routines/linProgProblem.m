function [f,A,b,Aeq,beq,intcon,uubLin,llbLin,param] = linProgProblem(param)
% linProgProblem  Build the MILP approximation used in mainMILP.
%
% Decision vector:
%   [dv_plus(3N); dv_minus(3N); z(1:disc-1)]
% where z are binary selectors for tangent half-planes.
%
% Input:
%   param   - struct containing linear maps, ellipse samples, dv bounds,
%             and variable dimensions
% Outputs:
%   f,A,b,Aeq,beq,intcon,uubLin,llbLin - MILP matrices/vectors for intlinprog
%   param   - returned unchanged for interface consistency

M = 1e4;
dum = [param.drrb2Nom -param.drrb2Nom];
for ind = 1:param.disc-1
    
    A(ind,:) = - param.xx2TargetNormal(:,ind)'*dum;
    b(ind,1) =   param.xx2TargetNormal(:,ind)'*(param.rrb2Nom-param.xx2Target(:,ind));

end

% Big-M activation of half-plane inequalities.
A = [A -M*eye(ind,ind)];

Aeq = [zeros(1,param.ndv*6) ones(1,ind)];
beq = ind-1;
A = [A;Aeq];
b = [b;beq];

% Penalize only continuous delta-v components.
f = zeros(1,param.nOptVar); f(1:end-ind) = 1000;

intcon = param.ndv*6+1:param.nOptVar;
uubLin =  [param.dvMax*ones(param.ndv*6,1); ones(ind,1)];
llbLin =  [zeros(param.ndv*6,1); zeros(ind,1)];

end

