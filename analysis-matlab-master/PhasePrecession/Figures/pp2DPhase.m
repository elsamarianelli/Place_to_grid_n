function [ret] = pp2DPhase(varargin)
% This is a quickly thrown together function that collects all the
% regressor values and phase values for either pdcd or pdmd and "split" and
% subset of cells and computes the regression for the group data, in
% addition it plots 2D phase and variance diagrams and linearised density
% and raster plots of phase versus pdcd or pdmd.
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Boilerplate for standardised input/output behaviour
if nargin == 0, ret = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

% choose correct regression structure.
if ~isempty(in.regress)
    regress = in.regress;
else
    regress = [in.inStruct.regress];
end

% get ind into regression parameters for the correct regressor/split
regressors = {'spks'; 'cfr'; 'ifr'; 'tim'; 'd'; 'pdmd'; 'pdcd'; 'tcn'};
whichReg = find(strcmp(regressors,in.regVar));

% initialise vars for the coming loop
nSpikes = zeros(numel(in.inStruct),1);
nCycles = zeros(numel(in.inStruct),1);
nPos = zeros(numel(in.inStruct),1);
intercept = zeros(numel(in.inStruct),1);
slope = zeros(numel(in.inStruct),1);

fnames = fieldnames(in.inStruct);

spk = vertcat(in.inStruct(:).spk);
pos = vertcat(in.inStruct(:).pos);
runs = vertcat(in.inStruct(:).runs);
theta = rmfield(in.inStruct,fnames([1:17,19:22]));
in.inStruct = rmfield(in.inStruct,fnames(17:22));

% calculate new pos and spk run labels if necessary
for n=1:numel(in.inStruct)
    if isnan(spk(n).time)
        continue
    end
    if ~strcmp(in.splitVar,'all')
        % Define vector of split variable labels set in ppRegression
        splitters = {'spike'; 'speed'; 'rate'; 'peripheral'; 'tortuosity'};
        
        % choose vector to be used to split data according splitVar
        whichSplit = strcmp(splitters,in.splitVar);
        
        % Splitting variables to choose from
        runIsHi = [runs.runIsHiSpike, runs.runIsHiSpeed, runs.runIsHiRate,...
            runs.runIsPeripheral, runs.runIsTortuous];
       
        % choose variable by which to split data in to upper/lower half by
        % number of spikes
        si = runIsHi(:,whichSplit);
        
        % relabel spk and pos run labels according labels of remaing spikes
        spk(n).runLabel = applyFilterToLabels( si==in.hilo, spk(n).runLabel );
        pos(n).runLabel = applyFilterToLabels( si==in.hilo, pos(n).runLabel );

    end
    
    % Relabel runLabels if we only want to use 1st or last spike in theta cycle
    if strcmp(in.whichSpikeInTheta,'first')
        spk(n).runLabel(~spk(n).isFirstInTheta) = false;
    elseif strcmp(in.whichSpikeInTheta,'last')
        spk(n).runLabel(~spk(n).isLastInTheta) = false;
    end
    
    goodPhase = ~isnan(theta(n).phase(spk(n).eegInd));
    
    slope(n) = regress(n).slope(whichReg);
    intercept(n) = regress(n).intercept(whichReg);
    nCycles(n,1) = max(spk(n).thetaCycleLabel(spk(n).runLabel>0 & goodPhase));
    nSpikes(n,1) = sum(spk(n).runLabel>0);
    nPos(n,1) = sum(pos(n).runLabel>0);
end

% get spk and pos, start and end inds for each cell to correctly accumulate
% regressors and phase across all cells
if strcmp(in.whichSpikeInTheta,'mean')
    N = sum(nCycles);
    spkStartInds = cumsum([1;nCycles(1:end-1)]);
    spkEndInds = cumsum(nCycles);
    
else
    N = sum(nSpikes);
    spkStartInds = cumsum([1;nSpikes(1:end-1)]);
    spkEndInds = cumsum(nSpikes);
    
end
M = sum(nPos);
posStartInds = cumsum([1;nPos(1:end-1)]);
posEndInds = cumsum(nPos);

% initialise values for the coming loop
reg = nan(N,1);
pha = nan(N,1);
reg2 = nan(M,1);
pha2 = nan(M,1);
r = nan(N,1);
phi = nan(N,1);

% this loop accumulates the relevant values across all cells.
for n=1:numel(in.inStruct)
    if isnan(spk(n).time)
        continue
    end
    
    % matrix of eeg indices reshaped by 5xlength(pos) in order to easily
    % index those eeg indices corresponding to each pos index.
    eegInds = 1:ceil(numel(theta(n).phase)/5)*5;
    eegInds = reshape(eegInds, [5,ceil(numel(theta(n).phase)/5)]);
    
    regs = {1; 1; 1; 1; 1; pos(n).d_meandir; pos(n).d_currentdir; 1};
    phase = theta(n).phase(spk(n).eegInd(spk(n).runLabel>0));
    goodPhase = ~isnan(phase);
    
    phaseRealignment = slope(n)*in.interceptAlign + intercept(n);

    if strcmp(in.whichSpikeInTheta,'mean')
        %get values of phase at spike times
        cycleLabels = spk(n).thetaCycleLabel(spk(n).runLabel>0);
        
        %get mean phase over theta cycles
        cycleComplexPhase = accumarray(cycleLabels(goodPhase),...
            exp(1i*double(phase(goodPhase))), [max(cycleLabels), 1], [], NaN);
        phase = angle(cycleComplexPhase); %circ_mean
        pha(spkStartInds(n):spkEndInds(n)) = mod(phase - phaseRealignment,2*pi);
        
        %get spk counts per theta cycle
        spkCountPerCycle = accumarray(cycleLabels(goodPhase), 1,...
            [max(cycleLabels),1], [], NaN);
        regressor = regs{whichReg}(spk(n).posInd(spk(n).runLabel>0));
        
        %get mean of regressor values over theta cycle where spikes fired
        regressor = accumarray(cycleLabels(goodPhase),...
            double(regressor(goodPhase)), [max(cycleLabels), 1], [], NaN)...
            ./spkCountPerCycle;
        reg(spkStartInds(n):spkEndInds(n)) = regressor; %Inds won't match up if there are nans in phase
        
        %
        rr = pos(n).r(spk(n).posInd(spk(n).runLabel>0));
        phiphi = pos(n).phi(spk(n).posInd(spk(n).runLabel>0));
        
        rr = accumarray(cycleLabels(goodPhase), double(rr(goodPhase)),...
            [max(cycleLabels), 1], [], NaN)./spkCountPerCycle;
        
        phiphi = accumarray(cycleLabels(goodPhase), double(phiphi(goodPhase)),...
            [max(cycleLabels), 1], [], NaN)./spkCountPerCycle;

        %get values of regressors at each pos within run
        r(spkStartInds(n):spkEndInds(n)) = rr;
        phi(spkStartInds(n):spkEndInds(n)) = phiphi;
    else
        %get values of phase at spike times
        phase = theta(n).phase(spk(n).eegInd(spk(n).runLabel>0));
        pha(spkStartInds(n):spkEndInds(n)) =  mod(phase - phaseRealignment,2*pi);
        %get values of regressors at spike times
        reg(spkStartInds(n):spkEndInds(n)) = regs{whichReg}(spk(n).posInd(spk(n).runLabel>0));
        
        %get values of regressors at each pos within run
        r(spkStartInds(n):spkEndInds(n)) = pos(n).r(spk(n).posInd(spk(n).runLabel>0));
        phi(spkStartInds(n):spkEndInds(n)) = pos(n).phi(spk(n).posInd(spk(n).runLabel>0));
    end
    
    %get values of phase at each pos within run
    eegRunLabels = eegInds(:,pos(n).runLabel>0);
    phase2 = theta(n).phase(eegRunLabels);
    phase2 = angle(nansum(exp(1i*phase2)));
    pha2((posStartInds(n)):posEndInds(n)) = mod(phase2 - phaseRealignment,2*pi);
    
    reg2(posStartInds(n):posEndInds(n)) = regs{whichReg}(pos(n).runLabel>0);
end
ret.pha3 = pha;
%discard NaNs
ind = isnan(reg) | isnan(pha) | isnan(r) | isnan(phi);
pha(ind) = [];
reg(ind) = [];
phi(ind) = [];
r(ind) = [];

% perform regression on grouped data
[regress_out.slope, regress_out.intercept] = circRegress(double(reg),double(pha));

theta = mod(abs(regress_out.slope)*double(reg), 2*pi);
[regress_out.cor, regress_out.p, regress_out.cor_boot,...
    regress_out.p_shuffled, regress_out.cor_ci] = ...
    circCorrTLinear(theta,double(pha),1000, 0.05, 0, true);

if in.plotFig
    % plot 2D phase and variance diagrams
    [x y] = pol2cart(phi,r);
    figure(length(get(0,'Children'))+1)
    in.saturation = 'even';
    spatialRotationalAverage(x,y,pha,[],in);
    title(['Phase map of average field, normalised to the unit circle, N = ',num2str(numel(in.inStruct))],'FontSize', 16);
    figure(length(get(0,'Children'))+1)
    in.saturation = 'variance';
    spatialRotationalAverage(x,y,pha,[],in);
    title('Phase variance map of average field, normalised to the unit circle','FontSize', 16);
    colorbar
    
    % and colour key if required
    if in.plotColorWheel
        figure(length(get(0,'Children'))+1)
        in.saturation = 'legend';
        spatialRotationalAverage(in);
        title('Colour wheel (Key to phase map)','FontSize', 16);
    end
    
    % plot phase vs regressor firing density
    reg(reg<=-1) = -1;
    reg(reg>=1) = 1;
    reg(reg==0) = eps;
    reg2(reg2<=-1) = -1;
    reg2(reg2>=1) = 1;
    reg2(reg2==0) = eps;
    
    ind = isnan(pha2) | isnan(reg2);
    pha2(ind) = [];
    reg2(ind) = [];
    
    ss = single([(mod(pha2(:)*(180/pi),360) + eps) ((reg2(:)+1)*180 + eps)]);
    ss_S = single([(mod(pha(:)*(180/pi),360) + eps) ((reg(:)+1)*180 + eps)]);
    si = accumarray(ceil(ss_S),1,[360 360]);
    ti = accumarray(ceil(ss),1,[360 360])/250;
    
    density_out = si./ti;
    density_out(isnan(density_out) | isinf(density_out))=0;
    
    [x y] = meshgrid(-1:(1/180):1, -pi:(pi/180):pi);
    sigmaX = 1/15;
    sigmaY = (1*pi)/15;
    kernel = exp(-0.5*((x/sigmaX).^2 + (y/sigmaY).^2));
    
    figure(length(get(0,'Children'))+1)
    hold off
    density_out = imfilter([density_out; density_out; density_out; density_out], kernel);
    density_out = density_out(size(si,1)+1:size(si,1)*3,:);
    imagesc(density_out); axis xy square
    title('Occupancy normalised phase x position rate map', 'FontSize', 16);
    set(gca, 'Box', 'Off', 'Xtick', [1 size(density_out,2)/2, size(density_out,2),...
        (3/2)*size(density_out,2), 2*size(density_out,2)],...
        'YTick', [1 size(density_out,1)/2, size(density_out,1),...
        (3/2)*size(density_out,1), 2*size(density_out,1)],...
        'XTickLabel', [-1 0 1], 'YTickLabel', {'0', 'pi', '2pi', '3pi', '4pi'},...
        'FontSize', 16, 'TickLength', [0.035 0]);
    xlabel('Normalised in-field position')
    ylabel('Theta Phase (radians)')
    
    % plot phase versus regressor raster plots with regression line
    figure(length(get(0,'Children'))+1)
    hold off
    plot([reg; reg], [pha; pha+2*pi],'k.', 'MarkerSize', 4);
    axis square; hold on;
    plot([-1 1], [-1*regress_out.slope + regress_out.intercept,...
        regress_out.slope + regress_out.intercept], 'r', 'LineWidth', 2);
    
    plot([-1 1], [-1*regress_out.slope + regress_out.intercept + ...
        2*pi, regress_out.slope + regress_out.intercept + 2*pi],...
        'r', 'LineWidth', 2);
    
    plot([-1 1], [-1*regress_out.slope + regress_out.intercept - ...
        2*pi, regress_out.slope + regress_out.intercept - 2*pi],...
        'r', 'LineWidth', 2);
    
    plot([-1 1], [-1*regress_out.slope + regress_out.intercept + ...
        4*pi, regress_out.slope + regress_out.intercept + 4*pi],...
        'r', 'LineWidth', 2);
    
    plot([-1 1], [-1*regress_out.slope + regress_out.intercept + ...
        6*pi, regress_out.slope + regress_out.intercept + 6*pi],...
        'r', 'LineWidth', 2);
    
    title(sprintf('Circular R = %4.3f | p-value < %4.3f',...
        regress_out.cor, regress_out.p_shuffled), 'FontSize', 16);
    set(gca, 'Xlim', [-1 1], 'YLim', [0 4*pi], 'Xtick', [-1 0 1],...
        'YTick', 0:pi:4*pi, 'XTickLabel', [-1 0 1], 'TickLength', [0.035 0],...
        'YTickLabel', {'0', 'pi', '2pi', '3pi', '4pi'}, 'LineWidth', 2, 'FontSize', 16);
    set(gca, 'Color', [204/255 1 1]);
    xlabel('Normalised in-field position')
    ylabel('Theta Phase (radians)')
end

ret.regress_out = regress_out;
% plot phase vs regressor averaged over each third of the field
% n = 3;
% reg = reg+1;
% reg(reg<=0) = eps;
% reg(reg>=2) = 2;
% 
% temp = zeros(n,3);
% 
% regress_out.cm = zeros(n,1);
% regress_out.cmci = zeros(n,2);
% for a = 1:n
%     t = (reg>(a-1)*(2/n) & reg<=(a*(2/n)));
%     regress_out.cm(a) = circ_mean(pha(t));
%     regress_out.cmci(a,1:2) = bootci(1000,{@circ_mean,pha(t)})';
%     temp(a) = numel(pha(t));
%     regress_out.pha(1:numel(pha(t)),a) = pha(t);
% end
% 
% if in.plotFig
%     cm = regress_out.cm;
%     cmci = regress_out.cmci;
%     figure(length(get(0,'Children'))+1)
%     X = -1+(2/(2*n)):(2/n):1;
%     e = errorbar(X(:),mod(cm-pi,2*pi),circ_dist(cm,cmci(:,1)),circ_dist(cm,cmci(:,2)));
%     h = findobj(gca,'Type','line');
%     set(h, 'LineWidth',3, 'Color', 'k')
%     axis square
%     box off
%     set(gca, 'Xlim', [-1 1], 'YLim', [-pi/16 2*pi+(pi/16)], 'Xtick', [-1 0 1],...
%         'YTick', 0:pi:4*pi, 'XTickLabel', [-1 0 1], 'TickLength', [0.035 0],...
%         'YTickLabel', {'0', 'pi', '2pi', '3pi', '4pi'}, 'LineWidth', 2, 'FontSize', 16);
%     xlabel('Normalised in-field position')
%     ylabel('Theta Phase (radians)')
%     
%     for er=1:length(e)
%         e2 = get(e(er), 'Children');
%         XData1 = get(e2(2), 'XData');
%         XData1(4:9:length(XData1)) = XData1(4:9:length(XData1)) - 0.04;
%         XData1(5:9:length(XData1)) = XData1(5:9:length(XData1)) + 0.04;
%         
%         XData1(7:9:length(XData1)) = XData1(7:9:length(XData1)) - 0.04;
%         XData1(8:9:length(XData1)) = XData1(8:9:length(XData1)) + 0.04;
%         
%         set(e2(2), 'XData', XData1);
%     end
% end
end

function in = DealWithInputs(varargin)
    defaults.inStruct = [];
    defaults.regress = [];
    
    defaults.regVar = 'pdcd';
    defaults.splitVar = 'all';
    defaults.whichSpikeInTheta = 'all';
    defaults.hilo = 0;
    defaults.interceptAlign = 0;
    defaults.plotColorWheel = false;
    defaults.plotFig = true;
    defaults.imageSize = 500;
    defaults.kernelSize = 50;
    defaults.kernelSigma = 10;
    defaults.saturation = 'even';

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