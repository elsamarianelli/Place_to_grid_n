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

%% Could fix alpha=1 (i.e. same amound of overlap as direction tuning)
% PS. not sure how good these will end up looking

alpha = 1; % amount of overlap between consecutive BVC's distal tuning 
s0 = 3; % baseline radial stdev
beta = 50; % distal sigma scaling

% s = x/beta +s0 is the radial stdev rule

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20

eqn1 = x1 == (x1/beta + s0)/alpha;
eqn2 = x2 == x1 + (x1/beta + x2/beta + 2*s0)/alpha;
eqn3 = x3 == x2 + (x2/beta + x3/beta + 2*s0)/alpha;
eqn4 = x4 == x3 + (x3/beta + x4/beta + 2*s0)/alpha;
eqn5 = x5 == x4 + (x4/beta + x5/beta + 2*s0)/alpha;
eqn6 = x6 == x5 + (x5/beta + x6/beta + 2*s0)/alpha;
eqn7 = x7 == x6 + (x6/beta + x7/beta + 2*s0)/alpha;
eqn8 = x8 == x7 + (x7/beta + x8/beta + 2*s0)/alpha;
eqn9 = x9 == x8 + (x8/beta + x9/beta + 2*s0)/alpha;
eqn10 = x10 == x9 + (x9/beta + x10/beta + 2*s0)/alpha;
eqn11 = x11 == x10 + (x10/beta + x11/beta + 2*s0)/alpha;
eqn12 = x12 == x11 + (x11/beta + x12/beta + 2*s0)/alpha;
eqn13 = x13 == x12 + (x12/beta + x13/beta + 2*s0)/alpha;
eqn14 = x14 == x13 + (x13/beta + x14/beta + 2*s0)/alpha;
eqn15 = x15 == x14 + (x14/beta + x15/beta + 2*s0)/alpha;
eqn16 = x16 == x15 + (x15/beta + x16/beta + 2*s0)/alpha;
eqn17 = x17 == x16 + (x16/beta + x17/beta + 2*s0)/alpha;
eqn18 = x18 == x17 + (x17/beta + x18/beta + 2*s0)/alpha;
eqn19 = x19 == x18 + (x18/beta + x19/beta + 2*s0)/alpha;
eqn20 = x20 == x19 + (x19/beta + x20/beta + 2*s0)/alpha;

sol = vpasolve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10,...
    eqn11, eqn12, eqn13, eqn14, eqn15, eqn16, eqn17, eqn18, eqn19, eqn20], ...
    [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, ...
    x12, x13, x14, x15, x16, x17, x18, x19, x20]);

sprintf('[%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f,%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f]', ...
    sol.x1,sol.x2,sol.x3,sol.x4,sol.x5,sol.x6,sol.x7,sol.x8,sol.x9,sol.x10, sol.x11, ...
    sol.x12,sol.x13,sol.x14,sol.x15,sol.x16,sol.x17,sol.x18,sol.x19,sol.x20)

%% Alternatively - could fix maximum distance (i.e. x20 = max_d)


s0 = 3; % baseline radial stdev
beta = 50; % distal sigma scaling
max_d = 177.3; % or whatever you want the maximum scaling to be
% Now we find alpha that gives the correct scaling up to the max distance
% s = x/beta +s0 is the radial stdev rule

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 alpha

eqn1 = x1 == (x1/beta + s0)/alpha;
eqn2 = x2 == x1 + (x1/beta + x2/beta + 2*s0)/alpha;
eqn3 = x3 == x2 + (x2/beta + x3/beta + 2*s0)/alpha;
eqn4 = x4 == x3 + (x3/beta + x4/beta + 2*s0)/alpha;
eqn5 = x5 == x4 + (x4/beta + x5/beta + 2*s0)/alpha;
eqn6 = x6 == x5 + (x5/beta + x6/beta + 2*s0)/alpha;
eqn7 = x7 == x6 + (x6/beta + x7/beta + 2*s0)/alpha;
eqn8 = x8 == x7 + (x7/beta + x8/beta + 2*s0)/alpha;
eqn9 = x9 == x8 + (x8/beta + x9/beta + 2*s0)/alpha;
eqn10 = x10 == x9 + (x9/beta + x10/beta + 2*s0)/alpha;
eqn11 = x11 == x10 + (x10/beta + x11/beta + 2*s0)/alpha;
eqn12 = x12 == x11 + (x11/beta + x12/beta + 2*s0)/alpha;
eqn13 = x13 == x12 + (x12/beta + x13/beta + 2*s0)/alpha;
eqn14 = x14 == x13 + (x13/beta + x14/beta + 2*s0)/alpha;
eqn15 = x15 == x14 + (x14/beta + x15/beta + 2*s0)/alpha;
eqn16 = x16 == x15 + (x15/beta + x16/beta + 2*s0)/alpha;
eqn17 = x17 == x16 + (x16/beta + x17/beta + 2*s0)/alpha;
eqn18 = x18 == x17 + (x17/beta + x18/beta + 2*s0)/alpha;
eqn19 = x19 == x18 + (x18/beta + x19/beta + 2*s0)/alpha;
eqn20 = max_d == x19 + (x19/beta + max_d/beta + 2*s0)/alpha;

sol = vpasolve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10,...
    eqn11, eqn12, eqn13, eqn14, eqn15, eqn16, eqn17, eqn18, eqn19, eqn20], ...
    [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, ...
    x12, x13, x14, x15, x16, x17, x18, x19, alpha]);

% This is a polynomial of order 20 and so there will be 20 possible solutions to this problem!
% Now we need to filter out the solution we're interested - should be able
% to do this by just looking for real, positive solutions (
sol_id = find( double(sol.x1) == real(double(sol.x1)) .* (real(double(sol.x1)) > 0) );

sprintf('[%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f]', ...
    sol.x1(sol_id),sol.x2(sol_id),sol.x3(sol_id),sol.x4(sol_id),sol.x5(sol_id),sol.x6(sol_id),sol.x7(sol_id),sol.x8(sol_id),sol.x9(sol_id),sol.x10(sol_id), sol.x11(sol_id), ...
    sol.x12(sol_id),sol.x13(sol_id),sol.x14(sol_id),sol.x15(sol_id),sol.x16(sol_id),sol.x17(sol_id),sol.x18(sol_id),sol.x19(sol_id), max_d)
