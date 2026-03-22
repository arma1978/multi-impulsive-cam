function [c,ceq,DC,DCeq] = nncon(xx,DAx,param)
% nncon  Nonlinear constraints and Jacobians for the NLP refinement stage.
%
% Inequalities: per-impulse norm constraints  ||dv_i|| <= dvMax
% Equality:     squared-Mahalanobis target constraint at encounter
%
% Inputs:
%   xx    - scaled optimisation vector [3*ndv x 1]
%   DAx   - DA maps loaded from runtime files
%   param - struct with scaling, bounds, reference trajectory, and maps
% Outputs:
%   c, ceq   - inequality and equality constraints for fmincon
%   DC, DCeq - analytical Jacobians of c and ceq (scaled variables)

xx = xx/param.dvScl;
%% extract data from param structure
dvMax = param.dvMax;
sqrMahaMin = param.sqrMahaMin;
rrvvNoMan = param.rrvvNoMan;

%% reshape maps for state evaluation
ndv = param.ndv;
dvv = reshape(xx, 3, ndv);
DV  = [zeros(size(dvv)); dvv; zeros(1,ndv)];
DAMaps = reshape(DAx(1:ndv*6),6, ndv);

DX = zeros(7, ndv);
dxf= zeros(6, 1);
rrvv = zeros(6,ndv);

%% evaluate the state after delta v applications
for ind =1:ndv
    dx0 = DX(:,ind)'+DV(:,ind)';
    for indj = 1:6
        dxf(indj,1)  = evalPoly(DAMaps(indj,ind).C, DAMaps(indj,ind).E, dx0);
    end
    DX(1:6,ind+1)= dxf-rrvvNoMan(:,ind);
    rrvv(:,ind) = dxf;
end

%% map the final state into conjuction and bplane to compute mahalanobis squared distance
DArrb2d = DAx(ndv*6+12+3:ndv*6+12+4);
DANb2drrb2d = DAx(ndv*6+12+5:ndv*6+12+6);

rrb2d(1,1) = evalPoly(DArrb2d(1).C, DArrb2d(1).E, DX(:,end)');
rrb2d(2,1) = evalPoly(DArrb2d(2).C, DArrb2d(2).E, DX(:,end)');
Nb2drrb2d(1,1) = evalPoly(DANb2drrb2d(1).C, DANb2drrb2d(2).E, DX(:,end)');
Nb2drrb2d(2,1) = evalPoly(DANb2drrb2d(2).C, DANb2drrb2d(2).E, DX(:,end)');
sqrMahalanobis = rrb2d'*Nb2drrb2d;
%% inequality and equality constraints
c1 = (sqrt(sum(DV.*DV,1))-dvMax)/dvMax;
c = c1';
ceq = (sqrMahaMin - sqrMahalanobis);

% Analytical Jacobians for fmincon (scaled variables)
DC = (cJacobian(dvv)/dvMax/param.dvScl)';
DCeq = -((xx'*param.Q'+ xx'*param.Q +param.L)/param.dvScl)';

end

function  DC = cJacobian(DDV)
% cJacobian  Jacobian of per-impulse norm constraints.

for ind = 1:size(DDV,2)
    DC(ind, (ind-1)*3+1:ind*3) = (DDV(:,ind))/norm(DDV(:,ind));
end

end