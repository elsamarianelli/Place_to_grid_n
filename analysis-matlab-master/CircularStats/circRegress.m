function [slope, intercept]=circRegress(x,t)
% Function to find approximation to circular-linear regression for phase
% precession.
% INPUTS:
%   x     - n-by-1 list of in-field positions (linear variable)
%   t     - n-by-1 list of phases
%         - neith can contain NaNs, must be paired (of equal length).
%   k     - no. of permutations to calculate p-value from randomisation  
%               test and confidence interval of correlation from bootstrap.
%   alpha - hypothesis test level e.g. 0.05 (5%) 0.01 (1%)
%   hyp   - hypothesis -1/0/1 [negatively correlated, correlated in either
%               direction, positively correlated]
%   conf  - true or false to calculate confidence intervals via jackknife
%           (highly optimised) or bootstrap - native matlab function
%
% Acknowledgments: 
% Fisher, N. I. 1996
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% NOTES:
% Circular-linear correlation needs to be taken out of this function. It
% should be run outside of it!

%% Example (Data from Fisher $3.3)
% x = [107 46 33 67 122 69 43 30 12 25 37 69 5 83 68 38 21 1 71 60 71 ...
%       71 57 53 38 70 7 48 7 21 27];
% t = [67 66 74 61 58 60 100 89 171 166 98 60 197 98 86 123 165 133 101 ...
%       105 71 84 75 98 83 71 74 91 38 200 56]*(pi/180);
% [slope, intercept]=circRegress(x,t)

%% Housekeeping
% Transform linear co-variate to lie between [-1 1].
mnx = mean(x);
x = x - mnx;
mxx = max(abs(x))+eps;
x = x/mxx;

% A not very robust check to keep angles in range [0 2*pi]
t = mod(t,2*pi);

% Check dimensional consistency of inputs
if any(size(x)~=size(t))
    error('x and t must be of the same size/dimensionality');
end

% Constrain maximum slope to give at most 720 deg phase precession over the
% field. 
max_slope = (4*pi)/(max(x)-min(x));

%% Performs slope optimisation using fminbnd and find intercept
[slope]=fminbnd(@(m) cost(m,x,t),...
    -1*max_slope, max_slope, optimset('TolFun',1e-3, 'TolX', 1e-2));

intercept = atan2(sum(sin(t - slope*x)), sum(cos(t - slope*x)));

intercept = intercept + (-slope*(mnx/mxx));
slope = slope/mxx;

%% Get Circular-Linear Association using Fisher / Jamalamakada
% theta = mod(abs(slope)*x, 2*pi);
% [cor p cor_b p_b ci] = circ_circ_corr(theta,t,k,alpha,hyp,conf);
% [cor p] = circ_corrcc(t, theta);
% cor_b = NaN;
% p_b = NaN;
% ci = [NaN NaN];
end

function [cost] = cost(m,x,t)
%% cost function to be minimised
%inputs: m=line slope, b=line intercept, x=data vector. cost is the
%resultant vector length.
alpha = t - m*x;
cost = -abs(sum(exp(1i*alpha)))/numel(alpha);

end