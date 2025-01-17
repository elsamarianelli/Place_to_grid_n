function varargout = plot_( x,y,varargin )
% Fast plotting function for showing long continuous signals.
% Inputs are exactly as for plot except...
% It currently ignores line style (: - .- etc.) and it doesn't support 
% extra input properties.
%
% x values must be in increasing order and close to being evenly spaced.
% 
% Example:
% x = 0:1e-3:1e3;
% plot_(x,sin(x),'k'); hold on; 
% plot_(x,cos(x),'r'); hold off;
%
% TODO: support more/all of the input forms that plot does
% TODO: do some checks on x to make sure it is suitable
% TODO: look into the fesibility of supporting tooltips.
%
% DM, Jan 2014.


if ishold %respect the current hold status
    xlim = get(gca,'XLim');
else
    cla;
    xlim = x([1 end]);
end

% Create the line object which will have it's x and y data set whenever we need to render the axes
extra= ParseExtraArgs(varargin{:});
hLine = line(x(1),y(1),extra{:});

% Store the full x and y data in the line's userdata
set(hLine,'UserData',struct('x',x,'y',y), 'DeleteFcn',@RemoveLineFromList);

% Add the line to the list of special lines within these axes
allLines = get(gca,'UserData');
if isempty(allLines)
    RegisterAxesCallback();  % see subfunction - this is what allows us to respond to axes changes
end
allLines(end+1) = hLine;
set(gca,'UserData',allLines);

% Draw the line for the current axes
RedrawLine(xlim,hLine,x,y);

if nargout == 1
    varargout{1} = hLine;
end
end

function extra = ParseExtraArgs(varargin)
extra = cell(0);
if mod(nargin,2) == 1
    linespec = varargin{1};
    colorInd = regexp(linespec,'[rgbcmyk]');
    
    extra(end+1:end+2) = {'Color',linespec(colorInd)};
end

end

% We assume x is sorted in increasing order and roughly evenly spaced.
function RedrawLine(xlims,hLine,x,y)
N = 1e3; %controls the resolution..ideally this should be roughly the number of pixels in the x-dimension

% Find the index that is just to the left of the screen
startInd = find(x>xlims(1),1,'first') -1;
if isempty(startInd)
    startInd = 1;
else
    startInd = max(startInd,1);
end

% And find the index that is just to the right of the screen
endInd = find(x>xlims(2),1,'first');
if isempty(endInd)
    endInd = numel(x);
end

% Calculate a stepSize so that there are roughly N divisions on the screen
stepSize = ceil((endInd-startInd+1)/N);
endInd = min(numel(x),endInd);

% Find the min and max y values in each of the N divisions
yBlock = y(startInd:endInd);
if numel(yBlock) < stepSize *N
    yBlock(end:stepSize*N) = yBlock(end);
end
yBlock = reshape(yBlock,stepSize,N);
[yBlockMax,maxInd] = max(yBlock,[],1);
[yBlockMin,minInd] = min(yBlock,[],1);

maxIsFirst = maxInd < minInd;

% Y data is interleaved as follows: y_left y_max y_min y_right, with the
% max and min swapped if they occured in the other order.
yVals([1 4],:) = yBlock([1 end],:);
yVals(2,maxIsFirst) = yBlockMax(maxIsFirst);
yVals(3,~maxIsFirst) = yBlockMax(~maxIsFirst);
yVals(3,maxIsFirst) = yBlockMin(maxIsFirst);
yVals(2,~maxIsFirst) = yBlockMin(~maxIsFirst);
yVals = yVals(:);

% Pick out, and duplciate, x values for each division
ind_a = (startInd:stepSize:endInd)';
ind_a(end+1:N) = ind_a(end);
ind_b = ind_a+stepSize-1;
ind_b(end) = min(ind_b(end),numel(x));
xVals = x([ind_a(:)'; ind_a(:)' ; ind_a(:)' ; ind_b(:)']);  
xVals = xVals(:);

% Done. Easy!
set(hLine,'XData',xVals,'YData',yVals);

end

function TriggerOnAllAxes(fg)
    ax = findall(fg,'type','axes');
    for ii=1:numel(ax)
        RedrawAllSeries(ax(ii));
    end
end

function RedrawAllSeries(ax)
theLines = get(ax,'UserData');
xlims = get(ax,'XLim');
for ii=1:numel(theLines)
    iiData = get(theLines(ii),'UserData');
    RedrawLine(xlims,theLines(ii),iiData.x,iiData.y);
end
end

function RemoveLineFromList(src,evt)
hLine = gcbo;
try
    ax = get(hLine,'Parent');
    allLines = get(ax,'UserData');
    allLines(allLines == hLine) = [];
    set(ax,'UserData',allLines);
    %TODO: I guess maybe we could check the length of allLines and
    %potneitally remove all the callbacks.
catch ex
end
end

function RegisterAxesCallback()
% This will call TriggerOnAllAxes whenever any axes are zoomed, panned, or
% resized.

set(zoom(gcf),'ActionPostCallback',@(src,evt) TriggerOnAllAxes(gcf));
set(pan(gcf),'ActionPostCallback',@(src,evt) TriggerOnAllAxes(gcf));
set(gcf,'ResizeFcn',@(src,evt) TriggerOnAllAxes(gcf));


end



