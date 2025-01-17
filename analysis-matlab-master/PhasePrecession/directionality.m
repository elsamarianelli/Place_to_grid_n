function [ret] = directionality(varargin)
%% Calculate extent of directional modulation of firing rate
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

%% Boilerplate for standardised input/output behaviour
if nargin == 0, ret = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

%% Compute Directionality and Skaggs Info
% Perhaps we should compare results to kldiv of firing rate versus
% distributive hypothesis instead of ml estimates.
spkPosInd = ceil(in.spikeTS*in.psr);
spkPosInd(spkPosInd > length(in.xy)) = size(in.xy,1);

% Convert data to bin units
xyBinUnits = toBinUnits([in.xy in.dir], [min(in.xy) 0; max(in.xy) 360], in.binsize);

% Create ratemap and smooth
xyBinInds = ceil(xyBinUnits);
dwell = accumarray(xyBinInds ,1)/in.psr;
spike = accumarray(xyBinInds(spkPosInd,:),1,size(dwell));

% Create ratemap
ml = maxLikelihoodRate(struct('dwell',dwell, 'spike',spike));
ret.converged = ml.converged;

if ~ml.converged
   spike = squeeze(sum(sum(spike,1),2));
   dwell = squeeze(sum(sum(dwell,1),2));
   nodwell = dwell==0;
   spike = imfilter(spike, in.dirKernel, 'circular');
   dwell = imfilter(dwell, in.dirKernel, 'circular');
   dirRate = spike./dwell;
   dirRate(nodwell) = interp1(dirRate(~nodwell), find(nodwell), 'linear', 'extrap');
else
    dirRate = imfilter(ml.ml_rate{2}, in.dirKernel, 'circular');
end
% % Get Skaggs Info (nb haven't done for pos ratemap)
% ti = squeeze(sum(sum(nd_dwell,1),2));
% [skInfo(1) skInfo(2)] = skaggs_info(ml.ml_rate{2}, ti);

% Get directionality
dirBinCentres = (360/60/2):(360/60):360;
% rayP = circ_rtest(dirBinCentres(:)*(pi/180), dirRate);
ret.divergence = kldiv(dirBinCentres(:)*(pi/180),ones(size(dirBinCentres(:)))/numel(dirBinCentres),dirRate(:)/sum(dirRate));

if in.PLOT_ON
    PlotForDebug(dirBinCentres(:),dirRate,ret.divergence)
end
end

function PlotForDebug(dirBinCentres,rate,kl)
% clf;
figure(1)
h = subplot(4,3,1);
polar((pi/180)*[dirBinCentres; dirBinCentres(1)],[rate; rate(1)])
hold on
% find all of the lines in the polar plot
L = findall(h,'type','line');                            % 1 - polar line
set(L([1 2 5 8]), 'LineWidth', 2, 'LineStyle', '-')      % 2:7 - radial grid-lines
set(L(10:2:end), 'LineWidth', 2, 'LineStyle', '-')       % 8 - outer concentric line
delete(L([3 4 6 7]))                                     % 9 - n will vary depending on max(rho) (keep alternate
delete(L(9:2:end))                                       % lines/text objects 9:2:end

% find all of the text objects in the polar plot
T = findall(gca,'type','text');                            % 1:3 - Don't know
set(T([4 5 10 11 16]), 'FontSize', 8)               % 4:15 - theta axis labels
set(T(18:2:end), 'FontSize', 8)                     % 16:end - rho axis labels
delete(T([1:2 6:9 12:15]))
delete(T(17:2:end))

title(sprintf('polar firing rate | %s', num2str(kl,2)), 'FontSize', 8)
set(h, 'FontSize', 8)
end

function in = DealWithInputs(varargin)

    defaults.xy = [];
    defaults.dir = [];
    defaults.psr = [];
    defaults.spikeTS = [];
    defaults.dirKernel = ones(5,1);
    defaults.binsize = [100/15 100/15 6];
    defaults.PLOT_ON = false;
    
    VERSION = 1;
    
    % Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0
       in = defaults;
       return;
    end
    if ischar(varargin{1})
        if mod(nargin,2) 
            error(['%s called with %d inputs, the first of which is a string. ' ...
                ' Either the first input should be a structure or there should be an even' ...
                ' number of inputs specifying structure field names and values.'],mfilename,nargin);
                
        end
        user_in = cell2struct(varargin(2:2:end)',varargin(1:2:end));
    else
        user_in = varargin{1};
        if nargin > 1 && VERSION ~= varargin{2}
            error(['%s called with version number %g, but code is version %g.\n' ...
                    'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
                    'the function without a version number.'],mfilename,varargin{2},VERSION);
        end
    end
    in = setstructfields(defaults,user_in);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end