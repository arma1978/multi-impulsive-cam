function param = ellipseTargets(param)
% ellipseTargets  Build B-plane risk geometry and target ellipse samples.
%
% Uses nominal encounter states/covariances from param to compute:
%   - baseline Mahalanobis distance and Pc estimates,
%   - target squared-Mahalanobis threshold,
%   - discretized ellipse points/normals used by convex solvers.

% param.rrvvsEncNoMan(1:3) = param.data{1,3:5}';
% param.rrvvsEncNoMan(4:6) = param.data{1,6:8}';
% param.rrvvdEncNoMan(1:3) = param.data{1,9+6:11+6}';
% param.rrvvdEncNoMan(4:6) = param.data{1,12+6:14+6}';

drr = param.rrvvsEncNoMan(1:3)-param.rrvvdEncNoMan(1:3);
dvv = param.rrvvsEncNoMan(4:6)-param.rrvvdEncNoMan(4:6);
Ps = param.PsRTN;
Pd = param.PdRTN;
%% RTN forward and backward transform
toRTNs = RTN(param.rrvvsEncNoMan);
toRTNd = RTN(param.rrvvdEncNoMan);
fromRTNd = toRTNd';
fromRTNs = toRTNs';

%% Covariances in a common inertial frame and Alfano Pc estimate
PsECI = fromRTNs*Ps*toRTNs;
PdECI = fromRTNd*Pd*toRTNd;
param.PcAlfano = PcAlfano( PsECI, PdECI,  ...
    drr, dvv, param.bodySize,100);
%% Project combined covariance to B-plane and evaluate baseline metrics
toBplane = Bplane(param.rrvvsEncNoMan(4:6), param.rrvvdEncNoMan(4:6));
Pb_3D = toBplane*(PsECI+PdECI)*toBplane';
Pb_2D(1,1) = Pb_3D(1,1); Pb_2D(1,2) = Pb_3D(1,3);
Pb_2D(2,1) = Pb_3D(3,1); Pb_2D(2,2) = Pb_3D(3,3);
Nb_2D = inv(Pb_2D);

dum = toBplane*drr;
drr2d = dum(1:2:3);
param.sqrMaha0 = drr2d'*Nb_2D*drr2d; 
param.maxPc0 = param.bodySize^2/exp(1)/param.sqrMaha0/det(Pb_2D)^0.5;
param.Pc0 = param.bodySize^2/2/det(Pb_2D)^0.5*exp(-param.sqrMaha0/2);

if isempty(param.Pc)
    param.sqrMahaMin = param.bodySize^2/exp(1)/param.PcMax/det(Pb_2D)^0.5;
else
    param.sqrMahaMin = -2*log(param.Pc*2*det(Pb_2D)^0.5/param.bodySize^2);
end

%% Ellipse of equal squared Mahalanobis distance
[V,D] = eig(Nb_2D);
[a,b] = sort(diag(D),'ascend');
D = diag(a);
V = V(:,b);

if isempty(param.minDist)
    a = sqrt(param.sqrMahaMin/D(1,1));
    b = sqrt(param.sqrMahaMin/D(2,2));
else
    a = param.minDist;
    b = a;
end

E = linspace(0,2*pi,param.disc);
xx2TargetEig = [a*cos(E);b*sin(E)];
tt = [-a*sin(E); b*cos(E)];
tt = tt./sqrt(sum(tt.*tt,1));
R = [0 1; -1 0];

%% Rotate ellipse samples back to the physical B-plane frame
% xx2Target = V*xx2TargetEig;
% xx2TargetNormal = V*R*tt;
param.xx2Target = V*xx2TargetEig;
param.xx2TargetNormal = V*R*tt;
param.ToBplaneCov = V';
param.semiMajor = a;
param.semiMinor = b;