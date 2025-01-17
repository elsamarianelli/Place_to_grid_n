function [regress] = ppRegression(theta,pos,spk,runs,regressionSet)
%% Calculate regression, correlation and p-values for multiple variables
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

% Redefine spk.runLabel as necessary for split data
newSpkRunLabel = spk.runLabel;
if ~strcmp(regressionSet.splitVar,'all')
    
    % Define vector of split variable labels set in ppRegression
    splitters = {'spike'; 'speed'; 'rate'; 'peripheral'; 'tortuosity'};
    
    % choose vector to be used to split data according splitVar
    whichSplit = strcmp(splitters,regressionSet.splitVar);
    runIsHi = [runs.runIsHiSpike, runs.runIsHiSpeed, runs.runIsHiRate,...
        runs.runIsPeripheral, runs.runIsTortuous];
    
    % choose variable by which to split data in to upper/lower half by
    % number of spikes
    si = runIsHi(:,whichSplit);
    
    % filter labels to reflect only those that are 'hilo'.
    newSpkRunLabel = applyFilterToLabels( si==regressioSet.hilo, newSpkRunLabel );
end

% Relabel runLabels if we only want to use 1st or last spike in theta cycle
spikeIsUsed = newSpkRunLabel>0;
if strcmp(regressionSet.whichSpikeInTheta,'first')
    spikeIsUsed(~spk.isFirstInTheta) = false;
elseif strcmp(regressionSet.whichSpikeInTheta,'last')
    spikeIsUsed(~spk.isLastInTheta) = false;
end

% Set regressor values
regressors = {spk.numWithinRun(spikeIsUsed);...
    pos.expectedRate_cumulative; pos.instantaneousFiringRate;...
    pos.timeSpentWithinRun; pos.d_cumulative; pos.d_meandir;...
    pos.d_currentdir;spk.thetaBatchLabelWithinRun(spikeIsUsed)};

% Get values of regressors at spike times
for i=2:7
    regressors{i} = regressors{i}(spk.posInd(spikeIsUsed));
end

% Get values of phase at spike times
phase = double(theta.phase(spk.eegInd(spikeIsUsed)));

% If we want to look at mean values within theta cycle:
if strcmp(regressionSet.whichSpikeInTheta, 'mean')

    goodPhase = ~isnan(phase);

    %get values of phase at spike times
    cycleLabels = spk.thetaCycleLabel(spikeIsUsed);
    
    %get mean phase over theta cycles
    cycleComplexPhase = accumarray(cycleLabels(goodPhase),...
        exp(1i*phase(goodPhase)), [max(cycleLabels), 1], [], NaN);
    phase = angle(cycleComplexPhase); %circ_mean
    
    %get spk counts per theta cycle
    spkCountPerCycle = accumarray(cycleLabels(goodPhase), 1,...
        [max(cycleLabels),1], [], NaN);

    for i=1:8
        %get mean of regressor values over theta cycle where spikes fired
        regressors{i} = accumarray(cycleLabels(goodPhase),...
            regressors{i}(goodPhase), [max(cycleLabels), 1], [], NaN)...
            ./spkCountPerCycle;
    end
end
goodPhase = ~isnan(phase);

% Initialise outputs
regress = struct('slope',nan(8,1), 'intercept',nan(8,1), 'cor',nan(8,1),...
    'p',nan(8,1), 'cor_boot',nan(8,1), 'p_shuffled',nan(8,1), 'ci',nan(8,2));

% Perform regression
for i=1:numel(regressors)
    goodRegressor = ~isnan(regressors{i});
    reg = regressors{i}(goodRegressor & goodPhase);
    pha = phase(goodRegressor & goodPhase);
    [regress.slope(i), regress.intercept(i)] = circRegress(reg, pha);
   
    mnx = mean(reg);
    reg = reg - mnx;
    mxx = max(abs(reg))+eps;
    reg = reg/mxx;
    
    theta = mod(abs(regress.slope(i))*reg, 2*pi);
    [regress.cor(i), regress.p(i), regress.cor_boot(i),...
        regress.p_shuffled(i), regress.ci(i,:)] = ...
        circCorrTLinear(theta,pha,regressionSet.k, regressionSet.alpha,...
        regressionSet.hyp, regressionSet.conf);
    
    % Option to calculate correlation using Jamalamdaka et al not fully
    % implemented:
    
    %     [regress.cor(i), regress.p(i), regress.p_shuffled(i)] = ...
    %         circCorrSRJ(theta,pha,regressionSet.k, regressionSet.hyp);

end

if regressionSet.PLOT_ON
    PlotForDebug(phase,regressors{7},regress.slope(7),...
        regress.intercept(7), regress.cor(7), regress.p_shuffled(7))
    PlotForDebug2(phase,regressors{6},regress.slope(6),...
        regress.intercept(6), regress.cor(6), regress.p_shuffled(6))
end

end

function PlotForDebug(phase,regressor,slope,intercept,cor,p_shuf)
figure(1)
subplot(4,3,10);
plot(regressor,phase,'k.'), axis square xy
hold on
plot(regressor,phase+2*pi,'k.'), axis square xy

% xl = get(gca,'XLim');
xl = [-1 1];

plot([-1 xl(2)], [-slope + intercept, slope + intercept], 'r', 'LineWidth', 2);
plot([-1 xl(2)], [-slope + intercept - 2*pi, slope + intercept - 2*pi],'r', 'LineWidth', 2);
plot([-1 xl(2)], [-slope + intercept + 2*pi, slope + intercept + 2*pi], 'r', 'LineWidth', 2);
plot([-1 xl(2)], [-slope + intercept + 4*pi, slope + intercept + 4*pi], 'r', 'LineWidth', 2);

set(gca, 'Xlim', [-1 1], 'YLim', [-pi 3*pi], 'Xtick', [-1 0 1],...
    'YTick', 0:pi:4*pi, 'XTickLabel', [-1 0 1],...
    'YTickLabel', {'0', 'pi', '2pi', '3pi', '4pi'}, 'FontSize', 6);
% set(gca, 'Xlim', [0 xl(2)], 'YLim', [-pi 3*pi], 'YTick', 0:pi:4*pi,...
%     'YTickLabel', {'0', 'pi', '2pi', '3pi', '4pi'}, 'FontSize', 6);
title(sprintf('phase versus normalised position \n S %s | I %s | C %s | P %s',...
    num2str(slope,2), num2str(intercept,2), num2str(cor,2), num2str(p_shuf,2)), 'FontSize', 8)
set(gca, 'FontSize', 8)
xlabel('pdcd');
ylabel('theta phase');
hold off;

end

function PlotForDebug2(phase,regressor,slope,intercept,cor,p_shuf)
figure(1)
subplot(4,3,11);
plot(regressor,phase,'k.'), axis square xy
hold on
plot(regressor,phase+2*pi,'k.'), axis square xy

% xl = get(gca,'XLim');
xl = [-1 1];

plot([-1 xl(2)], [-slope + intercept, slope + intercept], 'r', 'LineWidth', 2);
plot([-1 xl(2)], [-slope + intercept - 2*pi, slope + intercept - 2*pi],'r', 'LineWidth', 2);
plot([-1 xl(2)], [-slope + intercept + 2*pi, slope + intercept + 2*pi], 'r', 'LineWidth', 2);
plot([-1 xl(2)], [-slope + intercept + 4*pi, slope + intercept + 4*pi], 'r', 'LineWidth', 2);

set(gca, 'Xlim', [-1 1], 'YLim', [-pi 3*pi], 'Xtick', [-1 0 1],...
    'YTick', 0:pi:4*pi, 'XTickLabel', [-1 0 1],...
    'YTickLabel', {'0', 'pi', '2pi', '3pi', '4pi'}, 'FontSize', 6);
% set(gca, 'Xlim', [0 xl(2)], 'YLim', [-pi 3*pi], 'YTick', 0:pi:4*pi,...
%     'YTickLabel', {'0', 'pi', '2pi', '3pi', '4pi'}, 'FontSize', 6);
title(sprintf('phase versus normalised position \n S %s | I %s | C %s | P %s',...
    num2str(slope,2), num2str(intercept,2), num2str(cor,2), num2str(p_shuf,2)), 'FontSize', 8)
set(gca, 'FontSize', 8)
xlabel('pdmd');
ylabel('theta phase');
hold off;

end