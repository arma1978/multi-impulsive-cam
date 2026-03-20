function param = linearMaps(DAx, param)
% linearMaps  Extract first-order (linear) quantities from the DA flow maps.
%
%   Computes the state-transition matrix chain, the B-plane Jacobians,
%   and the quadratic/linear/constant terms of the squared Mahalanobis
%   distance as a function of the delta-v sequence.
%
%   Decision variable indexing: x = [ dv_1 (3); dv_2 (3); ... dv_N (3) ]
%   State transition chain: A(k) maps x -> deviation at node k
%   M = A(N): last-row block, maps x -> final state deviation  [6 x 3N]
%
%   Squared Mahalanobis distance (from B-plane geometry):
%     d^2 = dv'*Q*dv + 2*L*dv + c
%     where Q = M'*B'*C*M,  L = c1'*C*M + c2'*B*M,  c = c1'*c2
%
% Inputs:
%   DAx   - DA flow map cell array
%   param - struct with ndv
% Output:
%   param - updated with A, M, Q, L, c, drrb2Nom, rrb2Nom

ndv    = param.ndv;
DAMaps = reshape(DAx(1:ndv*6), 6, ndv);  % 6 x N

%% Build block state-transition matrix A and local STMs LL
%  A is (6N x 3N): A(k*6-5:k*6, :) maps the full dv vector to the
%  state deviation at node k.
LL = zeros(6, 6, ndv);
A  = zeros(ndv*6, ndv*3);
for ii = 1:ndv
    linear       = extractLinear(DAMaps(:,ii));  % linearise map at zero deviation
    LL(:,:,ii)   = linear(1:6, 1:6);             % 6x6 local STM (state wrt state)
    if ii == 1
        A(1:6, 1:3) = LL(:,4:6,ii);             % velocity columns of STM = impulse response
    else
        % Propagate previous columns forward, then add direct impulse at current node
        A((ii-1)*6+1:ii*6, :)              = LL(:,:,ii) * A((ii-2)*6+1:(ii-1)*6, :);
        A((ii-1)*6+1:ii*6, (ii-1)*3+1:ii*3) = LL(:,4:6,ii);
    end
end
param.A = sparse(A);
M = A((ii-1)*6+1:ii*6, :);  % [6 x 3N]  final-state deviation map

%% Extract B-plane linear maps from the DA maps
%  DArrb2d    : DA map of the 2D B-plane relative position  (c1 = constant, B = Jacobian)
%  DANb2drrb2d: DA map of Nb_2D * rrb2d                    (c2 = constant, C = Jacobian)
DArrb2d      = DAx(ndv*6+12+3 : ndv*6+12+4);
c1           = [DArrb2d(1).C(1); DArrb2d(2).C(1)];  % nominal B-plane position  [km]
linearDum    = extractLinear(DArrb2d);
B            = linearDum(:, 1:6);                    % [2 x 6]  d(rrb2d)/d(state)

DANb2drrb2d  = DAx(ndv*6+12+5 : ndv*6+12+6);
c2           = [DANb2drrb2d(1).C(1); DANb2drrb2d(2).C(1)];
linearDum    = extractLinear(DANb2drrb2d);
C            = linearDum(:, 1:6);                    % [2 x 6]  d(Nb2d*rrb2d)/d(state)

%% Quadratic form coefficients for squared Mahalanobis distance
Q = M' * B' * C * M;          % [3N x 3N]  quadratic coefficient matrix
L = c1'*C*M + c2'*B*M;        % [1  x 3N]  linear coefficient vector
c = c1' * c2;                 % scalar     constant term

% Store all quantities needed by linConvexProblem / linConvexProblemRefine
param.M       = M;        % [6 x 3N]  final-state deviation map
param.Q       = Q;        % quadratic coefficient   [km^-2]
param.L       = L;        % linear coefficient      [km^-1]
param.c       = c;        % constant term           [-]
param.drrb2Nom = B * M;   % [2 x 3N]  B-plane Jacobian w.r.t. dv vector
param.rrb2Nom  = c1;      % [2 x 1]   nominal B-plane position  [km]