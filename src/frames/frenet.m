function toFrenet = frenet(xx)
% frenet  Compute TNB (Frenet-like) frame transforms along a trajectory.
%
% Input:
%   xx - stacked states [6xN], rows: [r; v]
% Output:
%   toFrenet - transforms [3x3xN], rows [T; N; B]

rr = xx(1:3,:);
vv = xx(4:6,:);

for ind=1:size(xx,2)
    rr(:,ind) = rr(:,ind)/norm(rr(:,ind));
    tt(:,ind) = vv(:,ind)/norm(vv(:,ind));
    bb(:,ind) = cross(rr(:,ind),vv(:,ind));
    bb(:,ind) = bb(:,ind)/norm(bb(:,ind));
    nn(:,ind) = cross(bb(:,ind),tt(:,ind));
    nn(:,ind) = nn(:,ind)/norm(nn(:,ind));
    toFrenet(:,:,ind) = [tt(:,ind)'; nn(:,ind)'; bb(:,ind)'];
end
