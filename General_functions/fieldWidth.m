function fw = fieldWidth(x)
%FIELDWIDTH Summary of this function goes here
%   x is distance from boundary, outputs field width orthogonal to boundary
%   according to parameters fit by sander in his biggest environment

g = 3298.5;
h = 40;
w = 74.7;

% g = 3383.6;
% h = 36.5;
% w = 19;

fw = g*(1/h - h/(h^2 +x^2)) + w;

end
