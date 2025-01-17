function [F p dof] = watsonwilliams(s, sa)
%% Multisample testing for difference between mean angles
%  Watson-Williams test
%   Method followed is the same as in Chapter 26 Zar JH (1999), Biostatistical
%   Analysis Prentice Hall, ISBN: 013081542X, 9780130815422
% INPUTS:
%   s - list of angles
%   sa - list of integers corresponding to the groups the angles in s
%   belong t.
%
% OUTPUTS:
%   F - Critical value of F-statistic
%   dof - degrees of freedom for F-statistic
%   p - p-value for F-statistic
%
% Example 26.8 from Zar 1999, pp625-626.
% s1 = [94 65 45 52 38 47 73 82 90 40 87]*(pi/180);
% s2 = [77 70 61 45 50 35 48 65 36]*(pi/180);
% s = [s1(:);s2(:)];
% sa = [ones(size(s1(:)));2*ones(size(s2(:)))];
% [F p dof] = watsonwilliams(s, sa)
%
% Example 26.9 from Zar 1999, pp627-628.
% s1 = [135 145 125 140 165 170]*(pi/180);
% s2 = [150 130 175 190 180 220]*(pi/180);
% s3 = [140 165 185 180 125 175 140]*(pi/180);
% s = [s1(:);s2(:);s3(:)];
% sa = [ones(size(s1(:)));2*ones(size(s2(:)));3*ones(size(s3(:)))];
% [F p dof] = watsonwilliams(s, sa)
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


% Ensure angles s are between -pi and pi. (May not be necessary)
s = fixAngle(s);
% No. of sample means
k = max(sa);
N = numel(s);
R = N*circ_r(s);
% No. of observations N (in each (ideally) or in largest sample)
n = hist(sa,1:k);
anga = NaN*ones(k,max(n));
r = ones(k,1);
for i=1:k
    ind = sa==i;
    anga(i,1:n(i)) = s(ind);
    r(i) = n(i)*circ_r(s(ind));
end

rw = sum(r)/N;
checkAssumption(rw,N)
K = circ_r2kappa(rw);
K = 1+3/(8*K);    % correction factor
F = K*((((N-k)*(sum(r)-R)))/((k-1)*(N-sum(r))));
dof = [k-1 N-k];
p = 1-cdf('f',F,dof(1),dof(2));

function checkAssumption(rw,N)

if N > 10 && rw<.45
  warning('Test not applicable. Average resultant vector length < 0.45.') %#ok<WNTAG>
elseif N > 6 && rw<.5
  warning('Test not applicable. Average number of samples per population < 11 and average resultant vector length < 0.5.') %#ok<WNTAG>
elseif N >=5 && rw<.55
  warning('Test not applicable. Average number of samples per population < 7 and average resultant vector length < 0.55.') %#ok<WNTAG>
elseif N < 5
  warning('Test not applicable. Average number of samples per population < 5.') %#ok<WNTAG>
end
