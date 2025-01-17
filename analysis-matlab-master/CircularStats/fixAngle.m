function [x] = fixAngle(x)
% Ensure angles are between -pi and pi.
% x must be in radians
%
%%    Copyright (C) <2013>  <Daniel Manson> <d.manson@ucl.ac.uk> &
%%                          <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
%                           
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

x = mod(x+pi,2*pi)-pi;

end