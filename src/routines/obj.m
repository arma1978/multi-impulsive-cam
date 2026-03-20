function [fval,gval] = obj(xx,param)
% obj  NLP objective for minimum total delta-v.
%
% Inputs:
%   xx    - scaled optimisation vector [3N x 1]
%   param - scaling factors dvScl and objScl
% Outputs:
%   fval  - scaled objective value
%   gval  - scaled gradient

xx = xx/param.dvScl;

DDV = reshape(xx, 3, length(xx)/3);
fval = sum(sqrt(sum(DDV.*DDV,1)))*param.objScl;
gval = objgrad(DDV)*param.objScl/param.dvScl;

end

function gval = objgrad(DDV)
% objgrad  Gradient of sum(||dv_i||) wrt stacked components.

for ind = 1:size(DDV,2)
    gval((ind-1)*3+1:ind*3) = (DDV(:,ind))/norm(DDV(:,ind));
end

end