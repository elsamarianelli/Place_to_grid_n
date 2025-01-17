function [runs split] = runSplit(runs)
%% Split runs according to various behavioural and spiking criteria
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

% runs is a structure containing several vectors of length nRuns

% (a) Spikes fired
[runs.runIsHiSpike, split.spikeThresh, split.spikeSpikeCount]...
    = SplitCountsInHalfByValue(runs.spikeCount,runs.spikeCount);

% (b) Speed
[runs.runIsHiSpeed, split.speedThresh, split.speedSpikeCount]...
    = SplitCountsInHalfByValue(runs.meanSpeed,runs.spikeCount);

% (c) Firing rate.
[runs.runIsHiRate, split.rateThresh, split.rateSpikeCount]...
    = SplitCountsInHalfByValue(runs.rateInPosBins,runs.spikeCount);

% (d) Central vs Peripheral position in field.
[runs.runIsPeripheral split.peripheralThresh, split.peripheralSpikeCount]...
    = SplitCountsInHalfByValue(runs.centralPeripheral,runs.spikeCount);

% (e) Tortuosity
[runs.runIsTortuous, split.tortuosityThresh, split.tortuosityCount]...
    = SplitCountsInHalfByValue(runs.tortuosity,runs.spikeCount);
end


function [runIsHi, threshold, countTotals] = SplitCountsInHalfByValue(values,counts)

    [valuesSorted, sortingInds] = sort(values);
    countsSorted = counts(sortingInds);
    cumulativeCountsSorted = cumsum(countsSorted);
    cumulativeCountSortedFraction = cumulativeCountsSorted / cumulativeCountsSorted(end);
    [ignore, splitInd] = min(abs(0.5 - cumulativeCountSortedFraction));
    threshold = valuesSorted(splitInd);
    runIsHi = values > threshold;
    countTotals = [cumulativeCountsSorted(splitInd) ...
                   cumulativeCountsSorted(end)-cumulativeCountsSorted(splitInd)];
               
end