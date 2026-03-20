function tanomal = periods(sma,ecc,inc,argper,j2,req,mu)
% periods  Compute J2-perturbed anomalistic period.
%
% Inputs are classical elements and constants in consistent units.
% Output tanomal is the period between perigee passages.

pi2 = 2*pi;
% Keplerian period

tkepler = pi2 * sma * sqrt(sma / mu);


e = sin(inc) * sin(inc);

% perturbed mean motion


ar = (sma / req) * (sma / req);
sw = sin(argper);


% anomalistic period - time between perigee passages

tanomal = tkepler * (1.0 - 1.5 * j2 * (1.0 - 3.0 * e * sw * sw) ...
    / (ar * (1.0 - ecc) ^ 3));

