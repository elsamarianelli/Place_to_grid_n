function varargout = make_model_GC_module(xy, gridScale, gridOrient, dxy,...
            cellRate, drift_dth, drift_dxy)
% See readme.md in this folder, best viewed rendered with fonts and images
% on github.
%
% DM, May 2015.
%

posSampRate=50; % used for converting from ideal rate to expected spikes per pos samp


if ~exist('gridScale','var'), gridScale = 30; end;
if ~exist('gridOrient','var'), gridOrient = 10; end;
if ~exist('dxy','var'), dxy = [0 4; 6 2]; end;
n_cells = size(dxy,2);
if ~exist('cellRate','var'), cellRate = ones(n_cells,1); end;

true_xy = xy; % needed for plotting only
xy = xy + 1000; % add offset to keep xy positive (as required by lookup_ideal_rate            
if exist('drift_dth','var') && ~isempty(drift_dth)
    xy = cumrot(xy, -drift_dth); % this is complciated!!
end

if exist('drift_dxy','var') && ~isempty(drift_dxy)
    xy = xy + cumsum(drift_dxy, 1);
end

varargout = cell(n_cells,1);
for ii=1:n_cells
    rate = lookup_ideal_rate(xy, gridScale, gridOrient, dxy(:,ii)');
    spsPerPosSamp = poissrnd(rate/posSampRate*cellRate(ii));
    varargout{ii} = repeatInd(spsPerPosSamp);
end

PlotForDebug(true_xy, varargout{:})

end

function PlotForDebug(xy, varargin)
cla
plot(xy(:,1),xy(:,2),'k');
hold on

colors = 'bgrmyc';
for ii=1:length(varargin)
    plot(xy(varargin{ii},1),xy(varargin{ii},2),['s' colors(ii)],'MarkerFaceColor',colors(ii));
end
axis image; 
end

