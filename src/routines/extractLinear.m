function linear = extractLinear(DA)
% extractLinear  Extract first-order coefficients from DA polynomials.
%
% Input:
%   DA - DA polynomial array with fields E (exponents) and C (coefficients)
% Output:
%   linear - matrix of linear coefficients for each polynomial entry

for ind = 1:length(DA)    
    order = sum(DA(ind).E,2);
    lin = DA(ind).E(order==1,:);
    coe = DA(ind).C(order==1);
    linear(ind,:) = sum(coe.*lin,1);
end