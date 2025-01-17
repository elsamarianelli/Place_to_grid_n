function [h, p, obs] = permtest(a,b)
%% HELP
% 
% Checks the significance of the direction of the difference between values
% in two equal-length (paired) vectors (a and b).
%
% obs - mean difference of randomly sign-reversed values
% h   - hypothesis null or exceeded
% p   - obeserved p-value.

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

pm = zeros(10000,1);
for i=1:10000
    ind = randi([1 numel(a)],ceil(numel(a)/2),1);
    c = a-b;
    c(ind) = -c(ind);
    pm(i) = mean(c);
end
obs = mean(a-b);
h = obs<prctile(pm,5);
p = (1 + sum(pm<obs))/10000;
