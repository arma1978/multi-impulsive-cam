function toBplane = Bplane(vvs,vvd)
% Bplane  Build orthonormal B-plane basis from encounter velocities.
%
% Inputs:
%   vvs - spacecraft velocity at encounter [3x1]
%   vvd - debris velocity at encounter     [3x1]
% Output:
%   toBplane - DCM-like transform [3x3] with rows [xi; eta; zeta]

relVel = vvs-vvd;
eta = relVel/norm(relVel);
CrossVel = cross(vvd,vvs);
xi = CrossVel/norm(CrossVel);
zeta = cross(xi,eta);
zeta = zeta/norm(zeta);
    
toBplane = [xi'; eta'; zeta'];

end