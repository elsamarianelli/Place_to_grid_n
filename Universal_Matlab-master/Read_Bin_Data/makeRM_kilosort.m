function[inputVar,ratemaps, trialData] = makeRM_kilosort(inputVar, filesToConcat)
%% Generates ratemaps from kilosort/phy output
​
% Loads raw data, linearises for use on Radial Arm Maze and forms arm specific RMs. 
% Converts into a struct similar to the one formed from klustakwik data
​
% Input:
%   - inputVar.dataFold
%   - inputVar.dataFile
%   - in.pos
​
% NB I haven't added in the stability tests for kilosort data yet
%   stability tests - if any of these are true, split the raw data into
%   sections before creating RMs to look for stability (different types)
%   - in.alternateStability
%   - in.sessionStability
%   - in.alternateMinStability
​
​
%   - in.includeCentres (include centre of maze in linearisation?)
​
% Output:
%   -trialData: struct with spike data, pos data, pos_sample index
%   -ratemaps: struct containing ratemaps for each cell, seperated into the
%   arms of an eight arm radial-arm-maze
​
%% Provide some variables
​
inputVar.fileRoot = inputVar.dataFold;
inputVar.dataFold = [inputVar.fileRoot filesep 'PhyOutput'];
​
inputVar.spatialBinSize = 2; %(cm)
inputVar.speedThresh = 3; % running speed threshold (cm/s)
inputVar.rmSmooth = 2; % RM gaussian smoothing kernel (bins)
inputVar.filterCorners = true; % eliminate firing around corners?
inputVar.pxlPerCm = 2.74; % needs to be measured in environment. True for SJS's RAM
inputVar.binSizePix = inputVar.pxlPerCm*inputVar.spatialBinSize; %if bins are 2cm, pxlPerBin = 2.74*2
inputVar.vars = default_varsEA(inputVar.dataFold);
inputVar.tets2use = 1:16;
​
%% Import numpy files from phy output
% NB: you will need the npy-Matlab repo to do this
​
spikeTimes = readNPY([inputVar.dataFold filesep 'spike_times.npy']); % length nSpikes, gives sample number - divide by sample rate for s
spikeClusters = readNPY([inputVar.dataFold filesep 'spike_clusters.npy']);
clusterInfo = tdfread([inputVar.dataFold filesep 'cluster_info.tsv']);
​
%% Read in pos data, concatenate and post-process 
​
[posData, pos_samp] = concatPos_kilo(inputVar.fileRoot, inputVar.vars, filesToConcat, spikeTimes);
sleep1 = strcmp(posData.trialNames,'S1');
sleep2 = strcmp(posData.trialNames,'S2');
sleepSessions = find(sleep1+sleep2);
if length(sleepSessions)==2
    PosSamples_preSleep = 1:posData.trialFinalInd(sleepSessions(1)-1);
    PosSamples_postSleep = posData.trialFinalInd(sleepSessions(1))+1:posData.trialFinalInd(sleepSessions(2)-1);
    activeSessionsPos = [PosSamples_preSleep,PosSamples_postSleep];
else
    if sum(sleep1)==1
        PosSamples_preSleep = 1:posData.trialFinalInd(sleepSessions(1)-1);
        PosSamples_postSleep = posData.trialFinalInd(sleepSessions(1))+1:posData.trialFinalInd(end);
        activeSessionsPos = [PosSamples_preSleep,PosSamples_postSleep];
    elseif sum(sleep2)==1
        PosSamples_preSleep = 1:posData.trialFinalInd(sleepSessions(1)-1);
        activeSessionsPos = [PosSamples_preSleep];
    end
end
​
​
allOfPos = 1:length(posData.xy);
%% Force kilosort files into old klustakwik data structure
​
goodClusters = clusterInfo.id(clusterInfo.group(:,1)=='g'); %search for all clusters marked as good. CAREFUL not to use other cluster grouping in phy.
spikesToKeepMask = ismember(spikeClusters,goodClusters); % mask for those spikes that come from a 'good' cluster
​
% openEphys gives each tetrode a seperate depth value
for t = 1:max(clusterInfo.depth)
    tetrodeClusters = clusterInfo.id(clusterInfo.depth==t);
    trialData.tetrode(t).tetrodeMask = ismember(spikeClusters,tetrodeClusters);
    trialData.tetrode(t).cut = spikeClusters(trialData.tetrode(t).tetrodeMask&spikesToKeepMask);
    trialData.tetrode(t).timestamp = double(spikeTimes(trialData.tetrode(t).tetrodeMask&spikesToKeepMask))./30000; %divide by sample rate and convert to double
    trialData.tetrode(t).pos_sample = pos_samp(trialData.tetrode(t).tetrodeMask&spikesToKeepMask);
    trialData.tetrode(t).id = t;
end
​
trialData.clusterInfo = clusterInfo; % this has all sorts of helpful things like spike amplitude, firing rate, channel etc.
​
trialData.pos = posData;
clear posData;
​
%% Pre-process the data
% eliminate low running speeds and linearise data. Create masks for
% direction/session etc.
​
speedMask = trialData.pos.speed > inputVar.speedThresh; % mask to only include pos samples where the animal is running 
[armData, trialData] = lineariseArmData(trialData, inputVar.dataFold, allOfPos,activeSessionsPos, inputVar.includeCentres); % Split the pos_samples into individual arms on the RAM
xRange = findXLims(armData, trialData); % find min and max x value for all arms so that RMs will be the same size
​
[runDirection] = directionFilter(trialData); % split pos samples into in-bound and out-bound runs
[activeSessions] = activeSessionFilter(trialData); % remove sleep sessions (don't want to use them when making RMs)
​
​
%% Bin the POS data
%  Get the dwell times
trialData.lin.header = trialData.pos.header;
for i = 1:length(armData)
    outArmMask = (armData{i}.mask==1) & (runDirection.outMask ==1) & (activeSessions ==1);
    inArmMask = (armData{i}.mask==1) & (runDirection.inMask ==1) & (activeSessions ==1);
    outArmMask = find(outArmMask);
    inArmMask = find(inArmMask);
    ratemaps.arm{i}.outRun.posBinData = cb_bin_data('dwell_time','position', inputVar.binSizePix, trialData.lin, outArmMask, xRange);
    ratemaps.arm{i}.inRun.posBinData = cb_bin_data('dwell_time','position', inputVar.binSizePix, trialData.lin, inArmMask, xRange);
end
​
% Cycle through each arm, then each cell and make the RM
for i = 1:length(armData)
    cellCount = 1;
    for t = 1 : length(trialData.tetrode)   
​
        % Find all cells on that channel
        cellIDs = nonzeros(unique(trialData.tetrode(1,t).cut));
​
        % Then, if this channel is to be analysed and there are cells...
        if ~isempty(cellIDs) && ismember(trialData.tetrode(t).id,inputVar.tets2use)
​
            % find the posIndex for this arm in each direction
            armOutIdx = allOfPos(armData{i}.mask&runDirection.outMask&speedMask & activeSessions);      
            armInIdx = allOfPos(armData{i}.mask&runDirection.inMask&speedMask & activeSessions); 
            
            % Loop through those cells...
​
            for c = 1 : length(cellIDs)
                    
                % Record the tetrode and cluster number, all position samples
                ratemaps.cellN{cellCount,1}.cellInd = [t cellIDs(c)];
                ratemaps.cellN{cellCount,1}.posSample  = nonzeros(trialData.tetrode(1,t).pos_sample(trialData.tetrode(1,t).cut==cellIDs(c)));
                
                % Add the cell waveform duration to the ratemap struct
                 %ratemaps.cellN{cellCount,1}.waveformDur = trialData.tetrode(t).cell(c).waveformDur_microS;
                
​
                % Separate the position samples by arm and run direction
                spikeMaskOut = ismember(ratemaps.cellN{cellCount,1}.posSample, armOutIdx);
                spikeMaskIn = ismember(ratemaps.cellN{cellCount,1}.posSample, armInIdx);
                ratemaps.cellN{cellCount,1}.arm{i}.outSpikes = ratemaps.cellN{cellCount,1}.posSample(spikeMaskOut);
                ratemaps.cellN{cellCount,1}.arm{i}.inSpikes = ratemaps.cellN{cellCount,1}.posSample(spikeMaskIn);
​
                % Bin spikes
                ratemaps.cellN{cellCount,1}.arm{i}.outRun.binSpks = cb_bin_data('spikes','position',inputVar.binSizePix, armData{i}, ratemaps.cellN{cellCount,1}.arm{i}.outSpikes, xRange);
                ratemaps.cellN{cellCount,1}.arm{i}.inRun.binSpks = cb_bin_data('spikes','position',inputVar.binSizePix, armData{i}, ratemaps.cellN{cellCount,1}.arm{i}.inSpikes, xRange);
​
                if sum(ratemaps.cellN{cellCount,1}.arm{i}.outRun.binSpks) == 0
                    ratemaps.cellN{cellCount,1}.arm{i}.outRun.binSpks = ratemaps.cellN{cellCount,1}.arm{i}.outRun.binSpks';
                end
                
                if sum(ratemaps.cellN{cellCount,1}.arm{i}.inRun.binSpks) == 0
                    ratemaps.cellN{cellCount,1}.arm{i}.inRun.binSpks = ratemaps.cellN{cellCount,1}.arm{i}.inRun.binSpks';
                end
                
                % Generate smoothed rate maps
                ratemaps.cellN{cellCount,1}.arm{i}.outRun.ratemap = make_smooth_ratemap(ratemaps.arm{i}.outRun.posBinData, ratemaps.cellN{cellCount,1}.arm{i}.outRun.binSpks, inputVar.rmSmooth, 'gaus', 'norm');            
                %ratemaps.cellN{cellCount,1}.arm{i}.outRun.ratemap = ratemaps.cellN{cellCount,1}.arm{i}.outRun.ratemap(~isnan(ratemaps.cellN{cellCount,1}.arm{i}.outRun.ratemap));
                ratemaps.cellN{cellCount,1}.arm{i}.outRun.peakR = max(ratemaps.cellN{cellCount,1}.arm{i}.outRun.ratemap);
                
                ratemaps.cellN{cellCount,1}.arm{i}.inRun.ratemap = make_smooth_ratemap(ratemaps.arm{i}.inRun.posBinData, ratemaps.cellN{cellCount,1}.arm{i}.inRun.binSpks, inputVar.rmSmooth, 'gaus', 'norm');            
                %ratemaps.cellN{cellCount,1}.arm{i}.inRun.ratemap = ratemaps.cellN{cellCount,1}.arm{i}.inRun.ratemap(~isnan(ratemaps.cellN{cellCount,1}.arm{i}.inRun.ratemap));
                ratemaps.cellN{cellCount,1}.arm{i}.inRun.peakR = max(ratemaps.cellN{cellCount,1}.arm{i}.inRun.ratemap);
         
                cellCount = cellCount + 1;
                
            end
    
        end
    end
end