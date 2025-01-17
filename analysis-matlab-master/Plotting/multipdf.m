function multipdf(newpath, flnm)
% This function requires Adobe distiller be installed. It will create a PDF
% of figures, one page for each of the corresponding open matlab figures. 
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk> &
%%                          <Tom Hartley>
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

% Need to locate RunFileX.PS which will be in a location something like and
% copy it to the location where you want to create the PDF.
copyfile('D:\Program Files (x86)\Adobe\Acrobat 9.0\Acrobat\Xtras\RunFileX.PS',[newpath flnm '.ps']);

% Then open it for appending:
fmpps=fopen([newpath flnm '.ps'],'a');

%Then in a loop, print your figures to eps, and add lines to the ps file:
% set(get(0,'Children'), 'PaperOrientation', 'landscape', 'PaperPosition', [0.6345    0.6345   28.4084   19.7150],)
figNo = sort(get(0,'Children'));
for t=1:length(figNo)
	figure(figNo(t));
	outfile=[newpath num2str(t) '.eps'];
	%for the ps file you need forward slashes...
	outfile2=strrep(outfile,'\','/');
    orient portrait;
	print('-depsc', outfile);
	fprintf(fmpps,'(%s) prun\n',outfile2);
    set(figNo(t), 'Visible', 'off');
end

fclose(fmpps);

% Then print the ps file to pdf using acrobat distiller:
system(['"D:\Program Files (x86)\Adobe\Acrobat 9.0\Acrobat\AcroDist.exe"',...
    ' /F /J "D:\Program Files (x86)\Adobe\Acrobat 9.0\Acrobat\Settings\',...
    'Smallest File Size.joboptions" /Q /N ', '"', [newpath flnm], '.ps"']);

end