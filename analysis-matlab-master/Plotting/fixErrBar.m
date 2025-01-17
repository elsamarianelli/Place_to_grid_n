function fixErrBar(h,correction)
% Correct the short error bars you get with errorbar function. correction
% is the amount by which you want to increase the width of the errorbar in
% units of the x-axis.
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk> &
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

for er=1:length(h)
    h2 = get(h(er), 'Children');
    XData1 = get(h2(2), 'XData');
    XData1(4:9:length(XData1)) = XData1(4:9:length(XData1)) - correction;
    XData1(5:9:length(XData1)) = XData1(5:9:length(XData1)) + correction;
    
    XData1(7:9:length(XData1)) = XData1(7:9:length(XData1)) - correction;
    XData1(8:9:length(XData1)) = XData1(8:9:length(XData1)) + correction;
    
    set(h2(2), 'XData', XData1);
end