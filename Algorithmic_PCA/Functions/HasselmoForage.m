function [Position,Direction] = HasselmoForage(env,n_samples)
%HASSELMOFORAGE Summary of this function goes here
%   Detailed explanation goes here
% Initialise parameters for velocity and camera 
mu = 0; sigma = 0.2; b = 16;
v = 20; Dir = rand*2*pi; dt = 1/50;
% Create random velocity samples
RandomTurn = normrnd(mu,sigma,[1,n_samples]);
RandomVel = min(raylrnd(b,[1,n_samples]),49);
% Allocate memory for x, y components of position and velocity
Position = zeros(n_samples,2); start = datasample(find(env.L == 2),1);
Direction = zeros(n_samples,1);
[j,i] = ind2sub([env.dim_y, env.dim_x],start); Position(1,:) = [i j]; Direction(1) = Dir;
for step = 2:n_samples
    [dWall, aWall] = minDistAngle(env,round(Position(step-1,:)),Dir);
     v = RandomVel(step);
    %update speed and turn angle
    if(dWall < 7 && abs(aWall) < pi/2)
        angle = sign(aWall)*(pi/2-abs(aWall)) + RandomTurn(step);
        v = v - 0.5*max([0,(v-5)]); %slow down
    else
        angle = RandomTurn(step);
    end
    % move
    Position(step,:) = Position(step-1,:) + v*dt*[cos(-Dir),sin(-Dir)];
    if (env.L(round(Position(step,2)),round(Position(step,1))) ~= 2)
        Position(step,:) = Position(step-1,:);
    end
    % turn
    Dir = Dir + angle;
    Direction(step) = Dir;
end
 
figure
imagesc(~env.map');
colormap gray
hold on
plot(Position(:,2),Position(:,1),'LineWidth',2)
title(sprintf('Time = %f mins', n_samples/(60*50)))
 
axis tight
set(gcf,'color','w');
set(gca,'FontSize',16)
set(gcf,'Position',[100,100,800,800]);
end