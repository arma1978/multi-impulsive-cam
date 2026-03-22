function [ CollProb, varargout ] = PcAlfano( COV1, COV2, deltar, deltav, OBJ, varargin  )
% PcAlfano  Compute collision probability using the Alfano approximation.
%
% Syntax:
%   CollProb = PcAlfano(COV1, COV2, deltar, deltav, OBJ)
%   CollProb = PcAlfano(COV1, COV2, deltar, deltav, OBJ, N)
%
% Inputs:
%   COV1, COV2 - 6x6 state covariance matrices at TCA (position block used)
%   deltar     - relative position vector at TCA [3x1]
%   deltav     - relative velocity vector at TCA [3x1]
%   OBJ        - combined hard-body radius
%   N          - optional number of integration terms
%                (if omitted, selected automatically)
%
% Outputs: %   CollProb   - collision probability estimate


% Written by Alessandro Morselli, 24/10/2013


if nargin==6
    n = varargin{1};
else
    n = 0;
end

%% Select position covariance only
C = COV1(1:3, 1:3) + COV2(1:3, 1:3);

%% Define Bplane unit vectors
Xb = deltar/norm(deltar);

Yb = cross(deltar, deltav);
Yb = Yb/norm(Yb);

%% Map covariance matrix onto B-plane
R_XbYb = [ Xb'; Yb'];

C_B = R_XbYb*C*R_XbYb';

%% Map relative position onto B-plane
deltar_B = R_XbYb*deltar;

%% Obtain ellipse orientation in B-plane
[E, lambda] = eig(C_B);
lambda = diag(lambda);

[a, I] = max(sqrt(lambda));
e_a = E(:, I);
e_a = e_a/norm(e_a);

[b, I] = min(sqrt(lambda));
e_b = E(:, I);
e_b = e_b/norm(e_b);

% phi_B = acos(dot(e_a,Xb));

%% Rotate to covariance frame
R_ab = [e_a'; e_b'];
C_ab = R_ab*C_B*R_ab';

%% Map relative position on the B-plane
deltar_ab = R_ab*deltar_B;

xm = deltar_ab(1);
ym = deltar_ab(2);

sigmax = sqrt(C_ab(1,1));
sigmay = sqrt(C_ab(2,2));

%% Compute number of points for the integration
m = fix(5*OBJ/min([a,b,norm(deltar)]));

if n==0
    n = m;
end

if n<m
    fprintf('Warning: number of points for integration lower than suggested minimum\n');
    fprintf('n = %d\t, suggested at least m = %d\n',n,m);
end
    

CollProb = 0;
for i = 1:n
    
    aux1 = ( ym + 2*OBJ/n*sqrt((n-i)*i) )/(sigmay*sqrt(2));
    
    aux2 = (-ym + 2*OBJ/n*sqrt((n-i)*i) )/(sigmay*sqrt(2));
    
    aux3 = -(OBJ*(2*i-n)/n + xm )^2/(2*sigmax^2);
    
    CollProb = CollProb + (erf(aux1)+erf(aux2))*exp(aux3);
    
end

CollProb = OBJ*2/(sqrt(8*pi)*sigmax*n)*CollProb;

if nargout>1
    varargout{1} = C_B;
    if nargout>2
        varargout{2} = R_XbYb;
        if nargout>3
            varargout{3} = C_ab;
            if nargout>4
                varargout{4} = R_ab;
            end
        end
    end
end


end

