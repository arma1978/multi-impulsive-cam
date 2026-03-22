function p = evalPoly(C, E, x)
% evalPoly  Evaluate a DA/Taylor polynomial map at one or more points.
%
%   p = evalPoly(C, E, x)
%
%   Inputs:
%     C - coefficient column vector (K x 1), one entry per monomial
%     E - exponent matrix (K x n), one row of exponents per monomial
%     x - evaluation points (N x n), each row is one point
%
% Outputs: %     p - column vector of evaluated polynomial values (N x 1)
%
%   The polynomial is evaluated as:
%     p(i) = sum_k  C(k) * prod_j  x(i,j)^E(k,j)

p = zeros(size(x, 1), 1);

for i = 1:size(x, 1)
    if isempty(C)
        p(i) = 0;
    else
        % Replicate the point to align with all monomials, then raise to exponents
        xi   = ones(length(C), 1) * x(i, :);
        p(i) = sum(C .* prod(xi .^ E, 2));
    end
end
