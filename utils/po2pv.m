function [rr, vv] = po2pv(PO, mu)
% po2pv  Convert classical orbital elements to Cartesian position and velocity.
%
%   [rr, vv] = po2pv(PO, mu)
%
%   Inputs:
%     PO(1) - semi-major axis   a       [any consistent length unit]
%     PO(2) - eccentricity      e       [-]
%     PO(3) - inclination       i       [rad]
%     PO(4) - RAAN              Omega   [rad]
%     PO(5) - argument of perigee  omega [rad]
%     PO(6) - true anomaly      theta   [rad]
%     mu    - gravitational parameter   [length^3 / time^2]
%
%   Outputs:
%     rr    - position vector  (3x1)  [same length unit as a]
%     vv    - velocity vector  (3x1)  [length/time consistent with mu]

a = PO(1);
e = PO(2);
i = PO(3);
Om = PO(4);
om = PO(5);
theta = PO(6);

A = [cos(om+theta) -sin(om+theta) 0;
     sin(om+theta) cos(om+theta)  0;
                0             0   1];

% Handle negative inclination convention
if i < 0
    i = pi + i;
end
    
B = [1      0       0;
     0 cos(i)  -sin(i);
     0 sin(i)   cos(i)];
 
C =  [cos(Om) -sin(Om) 0;
      sin(Om) cos(Om)  0;
          0       0   1];
      
p = a*(1-e^2);

r = [1/(1+e*cos(theta))*p 0 0]';
v = sqrt(mu/p)*[e*sin(theta) 1+e*cos(theta) 0]';

rr = C*B*A*r;
vv = C*B*A*v;