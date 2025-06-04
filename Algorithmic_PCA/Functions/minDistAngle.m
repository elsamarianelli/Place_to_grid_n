function [dWall, aWall] = minDistAngle(env,pos,dir)
%MINDISTANGLE Returns distance and angle to the nearest wall.
%
%   [dWall, aWall] = minDistAngle(env, pos, dir)
%
%   Inputs:
%     env  - Environment struct (includes env.L, env.dwmap, env.map)
%     pos  - [x, y] position of the agent
%     dir  - Current heading direction (in radians)
%
%   Outputs:
%     dWall - Distance to the nearest wall (in pixels)
%     aWall - Relative angle (in radians) from heading to the wall
%
%   Notes:
%     aWall ≈ 0   → wall ahead
%     aWall < 0   → wall to the right
%     aWall > 0   → wall to the left

% figure; imagesc(env.dwmap); hold on; plot(pos(1), pos(2),'.', 'MarkerSize', 20); hold off

x = pos(1); y = pos(2);
dmap = zeros(size(env.L)); % EM bug fix
dmap(y,x) = 1;
dmap = bwdist(dmap);
dmap(env.L==1) = NaN;

%distance to wall from current x,y
dwpmap = dmap;
dwpmap(env.dwmap~=0) = NaN;
%  figure
%  imagesc(dwpmap);colorbar

%angle from current x,y
[cc,rr] = meshgrid(1:size(env.map,1),1:size(env.map,2));    
amap = atan2(rr-y,cc-x) + pi;
amap(isnan(dwpmap)) = NaN;

dWall = min(dwpmap,[],'all');
aWall = dir-amap(datasample(find(dwpmap==dWall),1));
%   figure


end

