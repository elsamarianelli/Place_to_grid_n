%% Code for Theo and Volker help reduce number of parameters in BVC model simulations

% Want 20 distance tunings and and 16 directional tunings
% This script was orignally written to find smooth transitions between
% consecutive distance tunings in the BVC space.
% Standard environment is 100x100, stretched is 100x175 (or 100x200)

% Imagining BVC tuning fields as residing on a dartboard/pinwheel, it would
% be quite reasonable to assume equal tuning overlap (i.e. number of sd)
% between adjacent cells on the pinwheel

% Directional tuning sigma of 0.2 radians seems to be pretty standard for
% generating fields from the model that look like BVCs. This is basically
% the same as pi/16, which for 16 angular tunings means adjacent cells  
% preferred angle distribution will overlap 1 sd from either mean

%% Alternatively - could fix maximum distance (i.e. x20 = max_d)


s0 = 8; % baseline radial stdev
beta = 12; % distal sigma scaling
max_d = 48; % or whatever you want the maximum scaling to be
% Now we find a that gives the correct scaling up to the max distance
% s = x/beta +s0 is the radial stdev rule

syms x1 x2 x3 x4 a 

eqn1 = x1 == (x1/beta + s0)/a;
eqn2 = x2 == x1 + (x1/beta + x2/beta + 2*s0)/a;
eqn3 = x3 == x2 + (x2/beta + x3/beta + 2*s0)/a;
eqn4 = x4 == x3 + (x3/beta + x4/beta + 2*s0)/a;
eqn5 = max_d == x4 + (x4/beta + max_d/beta + 2*s0)/a;

sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, ], ...
    [x1, x2, x3, x4, a]);

% This is a polynomial of order 20 and so there will be 20 possible solutions to this problem!
% Now we need to filter out the solution we're interested - should be able
% to do this by just looking for real, positive solutions (
sol_id = find( double(sol.x1) == real(double(sol.x1)) .* (real(double(sol.x1)) > 0) );

sprintf('[%.1f, %.1f, %.1f, %.1f, %.1f]', ...
    sol.x1(sol_id),sol.x2(sol_id),sol.x3(sol_id),sol.x4(sol_id),max_d)