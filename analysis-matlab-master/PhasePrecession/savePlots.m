function savePlots(SETTINGS,n,t)
%% Save "debug" plots to file
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

%open in matlab for appending:
fmpps=fopen([SETTINGS.outPath SETTINGS.outFile '.ps'],'a');
%set filename
figure(n);
fname=[SETTINGS.outPath num2str(t) '.eps'];
%for the ps file you need forward slashes...
fname2=strrep(fname,'\','/');
orient portrait;
%print to eps file
print(['-f',num2str(n)],SETTINGS.outputFigureFormat, fname);
%edit "master" ps file
fprintf(fmpps,'(%s) prun\n',fname2);
fclose(fmpps);

end