function [xxOnEllipse,xxOnEllipseNormal] = pointOnEllipse(xx0,param)
% pointOnEllipse  Nearest sampled point on target ellipse and outward normal.
%
% Inputs:
%   xx0   - query point in B-plane [2x1]
%   param - ellipse parameters (semi axes and frame)
% Outputs:
%   xxOnEllipse       - nearest sampled point [2x1]
%   xxOnEllipseNormal - corresponding unit normal [2x1]

N = 100000;
E = linspace(0,2*pi,N);
if isempty(param.minDist)
    a = param.semiMajor;
    b = param.semiMinor;
else
    a = param.minDist;
    b = a;
end

xx = param.ToBplaneCov'*[a*cos(E);b*sin(E)];

dxx = xx0-xx;
dx = sqrt(sum(dxx.*dxx,1));
[~,ind] = min(dx);
xxOnEllipse = xx(:,ind);

tt = [-a*sin(E(ind)); b*cos(E(ind))];
tt = tt/norm(tt);
R = [0 1; -1 0];
xxOnEllipseNormal = param.ToBplaneCov'*R*tt;

end
