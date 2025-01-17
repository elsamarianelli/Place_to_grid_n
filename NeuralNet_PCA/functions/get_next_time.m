function [x, y, direction] = get_next_time(x, y, In, direction)
%% trajectory generation used in dordek
%  do in funciton instead of out, doesnt have the boundary avoidance 
%  aspect is the only problem so try incorporate this later

angular_velocity = .2; % try changing? (org = .2)
direction = mod(direction + angular_velocity*randn,2*pi);

% next x after modulo (for spatial continuity)
next_x = mod(x+angular_velocity*cos(direction),In.x);
next_y = mod(y+angular_velocity*sin(direction),In.y);

% getting new location
x = floor(next_x);
y = floor(next_y);

end