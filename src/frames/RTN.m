function toRTN = RTN(xx)
% RTN  Compute RTN frame transforms along a trajectory.
%
% Input:
%   xx - stacked states [6xN], rows: [r; v]
% Output:
%   toRTN - transforms [3x3xN], rows [R; T; N]

rr = xx(1:3,:);
vv = xx(4:6,:);

for ind=1:size(xx,2)
    rr(:,ind) = rr(:,ind)/norm(rr(:,ind));
    vv(:,ind) = vv(:,ind)/norm(vv(:,ind));
    hh(:,ind) = cross(rr(:,ind),vv(:,ind));
    hh(:,ind) = hh(:,ind)/norm(hh(:,ind));
    tt(:,ind) = cross(hh(:,ind),rr(:,ind));
    tt(:,ind) = tt(:,ind)/norm(tt(:,ind));
    toRTN(:,:,ind) = [rr(:,ind)'; tt(:,ind)'; hh(:,ind)'];
end
