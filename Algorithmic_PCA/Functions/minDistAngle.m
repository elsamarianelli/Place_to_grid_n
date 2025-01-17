function [dWall, aWall] = minDistAngle(env,pos,dir)
%MINDISTANGLE Summary of this function goes here
%   Detailed explanation goes here
x = pos(1); y = pos(2);
dmap = zeros(env.dim_y,env.dim_x);
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

