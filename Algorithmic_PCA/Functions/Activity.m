function [ phi ] = Activity( x, y, cells, env, max_x,min_x,max_y,min_y)
%ACTIVITY Summary of this function goes here
%   Detailed explanation goes here

%distance from current x,y
dmap = zeros(env.dim_y,env.dim_x);
dmap(y,x) = 1;
dmap = bwdist(dmap);
dmap(env.L==1) = NaN;
%   figure
%   imagesc(dmap);colorbar

%distance to wall from current x,y
dwpmap = dmap;
dwpmap(env.dwmap>0) = NaN;
imap = nan(env.dim_y,env.dim_x);
imap(min_y:max_y, min_x:max_x) = 1; 
dwpmap = dwpmap.*imap;
%   figure
%   imagesc(dwpmap);colorbar

%angle from current x,y
[cc,rr] = meshgrid(1:size(env.map',1),1:size(env.map',2));    
amap = wrapTo2Pi(atan2(rr-y,cc-x));
amap(isnan(dwpmap)) = NaN;
%   figure
%   imagesc(amap);colorbar


phi = FiRate(dwpmap, amap, cells);
end

