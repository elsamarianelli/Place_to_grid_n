function [SETTINGS] = ppSettings(outFile)
%% Function to produce all neccesary constants for ppMain and functions
%% called therein [as input to ppMain]
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

close all hidden

%% Housekeeping for ppMain itself
% path and unique filename for output data.
SETTINGS.outPath = 'E:\Science\Matlab Programming\Phase Precession\ppOut\';
SETTINGS.outFile = [outFile,' ',datestr(now,'dd-mmm-HHMM')];

% This must be EPS if you are going to use distiller to create a single
% large PDF, otherwise it can be any figure format that matlab supports,
% one figure will be created per cell in SETTINGS.outPath and will be name
% by the numeral of the loop in which it was created.
SETTINGS.outputFigureFormat = '-depsc2';
% You can only do this if you have the full version of Adobe Acrobat
SETTINGS.outputCollatedPDF = true;
% make a copy the "master" eps file in the location you want to edit
if SETTINGS.outputCollatedPDF
    copyfile('D:\Program Files (x86)\Adobe\Acrobat 9.0\Acrobat\Xtras\RunFileX.PS',...
        [SETTINGS.outPath SETTINGS.outFile '.ps']);
end

% Define theta band - used throughout
SETTINGS.cutOffFreq = [6 12];
% Define thresholds on which to ignore cells
SETTINGS.gridnessThreshold = 0;
SETTINGS.intrThetaModulationThreshold = 5;
SETTINGS.LFPThetaModulationThreshold = 10;

% Global variable to set all debug plotting on or off. To choose individual
% plots go to SETTINGS structure for the relevant function below.
debugPlot = true;

%% Provide inputs required for load function
% Create 6-12Hz bandpass filter using 251-tap Blackman window. Numerators
% are real band pass frequencies, denominator is the Nyquist frequency.
SETTINGS.PPLOAD.EEGFilter = fir1(250,[SETTINGS.cutOffFreq(1)/125 SETTINGS.cutOffFreq(2)/125],blackman(251));

%% Provide inputs required for gcProps
% Define input structures for functions called within gcProps used to
%  calculate cell properties
SETTINGS.GCPROPS.xy = [];                                                   % Set within gcProps - xy position
SETTINGS.GCPROPS.psr = [];                                                  % Set within gcProps - position sampling rate
SETTINGS.GCPROPS.spikeTS = [];                                              % Set within gcProps - spike timestamps
SETTINGS.GCPROPS.sac = [];                                                  % Set within gcProps - spatial autocorrelation

SETTINGS.GCPROPS.binsize = [2 2];                                           % Binsize in cm

% autoCorr2D input
SETTINGS.GCPROPS.ac2D_in.x = [];                                            % Set within gcProps
SETTINGS.GCPROPS.ac2D_in.nodwell = [];                                      % Set within gcProps

SETTINGS.GCPROPS.ac2D_in.tol = 10^-10;                                      % tolerance to account for tiny numerical errors using fft method
SETTINGS.GCPROPS.ac2D_in.PLOT_ON = debugPlot;                               % Plot debug figures

% Autocorrelation
SETTINGS.GCPROPS.acProps_in.FIND_CENTROID = false;              
SETTINGS.GCPROPS.acProps_in.GET_PERIM_FIELDS = false;
SETTINGS.GCPROPS.acProps_in.GET_PERIM_GRIDNESS_MASK = false;
SETTINGS.GCPROPS.acProps_in.GET_ORIENTATION = true;                         % These are the only two properties we are after
SETTINGS.GCPROPS.acProps_in.GET_MEAN_R_AND_SCALE = true;                    % These are the only two properties we are after
SETTINGS.GCPROPS.acProps_in.FULL_GRIDNESS_RANGE = false;
SETTINGS.GCPROPS.acProps_in.FIELD_EXTENT_METHOD = 2;
SETTINGS.GCPROPS.acProps_in.PLOT_ON = debugPlot;                            % Plot debug figures

%% Provide inputs required for thetaMod
% Minimum Speed
SETTINGS.THETAMOD.eeg = [];                                                 % Set within thetaMod - eeg signal
SETTINGS.THETAMOD.spikeTS = [];                                             % Set within thetaMod - spike timestamps
SETTINGS.THETAMOD.spd = [];                                                 % Set within thetaMod - speed of movement
SETTINGS.THETAMOD.psr = [];                                                 % Set within thetaMod - position sampling rate
SETTINGS.THETAMOD.esr = [];                                                 % Set within thetaMod - eeg sampling rate

SETTINGS.THETAMOD.minSpeed = 5;                                             % Minimum speed below which to ignore singal (in cm/s)
SETTINGS.THETAMOD.minDuration = 0.5;                                        % Shortest allowable length of eeg signal in seconds (after filtering for speed)

% EEG power spectrum
SETTINGS.THETAMOD.eeg_in.pos2use = [];                                      % Set within eegPowerSpectra - valid position samples to use
SETTINGS.THETAMOD.eeg_in.eeg = [];                                          % Set within eegPowerSpectra - eeg signal
SETTINGS.THETAMOD.eeg_in.sampFreq = [];                                     % defined by header data and passed into thetamod

SETTINGS.THETAMOD.eeg_in.thetaRange = SETTINGS.cutOffFreq;                  % theta band
SETTINGS.THETAMOD.eeg_in.padToPow2 = NaN;                                   % defined by eeg signal
SETTINGS.THETAMOD.eeg_in.smthKernelWidth = 2;
SETTINGS.THETAMOD.eeg_in.smthKernelSigma = 0.1875;
SETTINGS.THETAMOD.eeg_in.s2nWdth = 2;
SETTINGS.THETAMOD.eeg_in.maxFreq = 125;
SETTINGS.THETAMOD.eeg_in.EEG_samples_per_position = 5;
SETTINGS.THETAMOD.eeg_in.xmax = 25;
SETTINGS.THETAMOD.eeg_in.ymax = [];
SETTINGS.THETAMOD.eeg_in.PLOT_ON = debugPlot;

% Spike autocorrelation power spectrum
SETTINGS.THETAMOD.ac_in.spikeTimes = [];                                    % Set within eegPowerSpectra - valid position samples to use
SETTINGS.THETAMOD.ac_in.posMask = [];                                       % Set within eegPowerSpectra - valid position samples to use
SETTINGS.THETAMOD.ac_in.thetaRange = SETTINGS.cutOffFreq;                   % theta band
SETTINGS.THETAMOD.ac_in.padToPow2 = 16;                                     % zero-pad spike intensity signal
SETTINGS.THETAMOD.ac_in.posSampFreq = 50;                                   % set by 
SETTINGS.THETAMOD.ac_in.acBinSize = 0.002;                                  % binsize in seconds for autocorrelation 
SETTINGS.THETAMOD.ac_in.acWindow = 0.5;                                     % maximum size in seconds of autocorrelation. Should be the same as SETTINGS.THETAMOD.minDuration
SETTINGS.THETAMOD.ac_in.smthKernelWidth = 2;                                % size in Hz of gaussian smoothing kernel
SETTINGS.THETAMOD.ac_in.smthKernelSigma = 0.1875;                           % standard deviation of gaussian smoothing kernel
SETTINGS.THETAMOD.ac_in.s2nWdth = 2;                                        % width around the peak used to calculate signal to noise ratio of signal in the theta band
SETTINGS.THETAMOD.ac_in.maxFreq = 125;                                      % frequency limit used of signal used as noise in s2n. Some people use 25 and 50.
SETTINGS.THETAMOD.eeg_in.xmax = 25;
SETTINGS.THETAMOD.eeg_in.ymax = [];
SETTINGS.THETAMOD.ac_in.PLOT_ON = debugPlot;                                % Plot figures for debugging

%% Construct settings required for directionality
SETTINGS.DIRECTIONALITY.xy = [];                                            % Set within directionaliy/ppMain - position data
SETTINGS.DIRECTIONALITY.dir = [];                                           % Set within directionaliy/ppMain - direction data
SETTINGS.DIRECTIONALITY.psr = [];                                           % Set within directionaliy/ppMain - position sample rate
SETTINGS.DIRECTIONALITY.spikeTS = [];                                       % Set within directionaliy/ppMain - spike timestamps
SETTINGS.DIRECTIONALITY.binsize = [100/15 100/15 6];                        % binsizes for position x in cm by position y in cm by direction in degress for pxd analysis of directionality

% smoothing kernel for directionality
dir_kernel_length = 180; %degrees
smthKernelSigma = 15; %degrees
binsPerDegree = 1/6;
kernelLengthInBins = round(dir_kernel_length*binsPerDegree);
kernelSigmaInBins = smthKernelSigma*binsPerDegree;
SETTINGS.DIRECTIONALITY.dirKernel = fspecial('gaussian',[kernelLengthInBins 1],kernelSigmaInBins);
SETTINGS.DIRECTIONALITY.PLOT_ON = debugPlot;

%% Construct settings required for getThetaProps
SETTINGS.GETTHETAPROPS.spikeTS = [];                                        % Set within getThetaProps - spike timestamps
SETTINGS.GETTHETAPROPS.phase = [];                                          % Set within getThetaProps - eeg phase
SETTINGS.GETTHETAPROPS.filteredEEG = [];                                    % Set within getThetaProps - filtered eeg signal
SETTINGS.GETTHETAPROPS.esr = [];

SETTINGS.GETTHETAPROPS.MEAN_POWER_PERCENTILE_THRESHOLD = 0;                 % percentile power below which to reject theta cycles
SETTINGS.GETTHETAPROPS.ALLOWED_THETA_CYCLE_LENGTH = [20 42];                % corresponds to SETTINGS.cutOffFreq in eeg sampling.
SETTINGS.GETTHETAPROPS.CANONICAL_MINIMUM_SPIKING_PHASE = pi;                % 

%% Construct settings required for partitionFields
SETTINGS.PARTITIONFIELDS.xy = [];                                           % Set within partitionFields - position data
SETTINGS.PARTITIONFIELDS.psr = [];                                          % Set within partitionFields - defined by header data and passed down
SETTINGS.PARTITIONFIELDS.spikeTS = [];                                      % Set within partitionFields - spike timestamps
SETTINGS.PARTITIONFIELDS.type = [];                                         % Set within partitionFields - cell type eg 'g' for grid, 'p' for place.

SETTINGS.PARTITIONFIELDS.binsize = [0.5 0.5];                               % binsize in cm
SETTINGS.PARTITIONFIELDS.FIELD_THRESHOLD_FRAC = 0.35;                       % this is a watershed input used to as a cutoff to define the field perimeter
SETTINGS.PARTITIONFIELDS.AREA_THRESHOLD_FRAC = 0.75;                        % used to discard fields that are too small and on the edge
% Define smoothing kernel
partition_kernel_length = 50; %centimetres
smthKernelSigma = 7.5; %centimetres
binsPerCm = 2;
kernelLengthInBins = round(partition_kernel_length*binsPerCm);
kernelSigmaInBins = smthKernelSigma*binsPerCm;
SETTINGS.PARTITIONFIELDS.smoothingKernelForPartiton = fspecial('gaussian',...
    [kernelLengthInBins kernelLengthInBins],kernelSigmaInBins);
SETTINGS.PARTITIONFIELDS.PLOT_ON = debugPlot;

%% Construct settings required for getPosProps        
SETTINGS.GETPOSPROPS.xy = [];                                               % Set within getPosProps - positions data in cm
SETTINGS.GETPOSPROPS.dir = [];                                              % Set within getPosProps - direction data
SETTINGS.GETPOSPROPS.spd = [];                                              % Set within getPosProps - speed data
SETTINGS.GETPOSPROPS.spikeTS = [];                                          % Set within getPosProps - spike time stamps 
SETTINGS.GETPOSPROPS.psr = [];                                              % Set within getPosProps - position sampling rate
SETTINGS.GETPOSPROPS.SPATIAL_LOW_PASS_CUTOFF = 3;

SETTINGS.GETPOSPROPS.peaksXY = [];                                          % Set within getPosProps - peak positions
SETTINGS.GETPOSPROPS.rmFieldLabel = [];                                     % Set within getPosProps - field labels

SETTINGS.GETPOSPROPS.binsize = [0.5 0.5];                                   % binsize in cm
SETTINGS.GETPOSPROPS.RUN_SMOOTH_WINDOW_FRACTION = 1/3;                      % IFR smoothing constant
SETTINGS.GETPOSPROPS.RUN_MINIMUM_DURATION = 2;                              % in pos bins.
SETTINGS.GETPOSPROPS.RUN_MINIMUM_SPEED = 2.5;                               % in cm/s
SETTINGS.GETPOSPROPS.smoothingWidthInPosBins = 15;                          % to smooth speed by

% smoothing kernel for instantaneous firing rate
ifr_kernel_length = 1; %seconds
smthKernelSigma = 0.5; %seconds
binsPerSec = 50;
kernelLengthInBins = round(ifr_kernel_length*binsPerSec);
kernelSigma = smthKernelSigma*binsPerSec;
SETTINGS.GETPOSPROPS.ifr_kernel = fspecial('gaussian',[kernelLengthInBins 1],kernelSigma);
SETTINGS.GETPOSPROPS.PLOT_ON = debugPlot;

%% Construct settings required for ppRegression
SETTINGS.PPREGRESSION.k = 1000;                                             % no of permutations for calculation of shuffle
SETTINGS.PPREGRESSION.alpha = 0.05;                                         % hypothesis test level
SETTINGS.PPREGRESSION.hyp = 0;                                              % single or two tailed test, -1, 0 (two-tailed), 1
SETTINGS.PPREGRESSION.conf = true;                                          % calculate confidence interval of correlation and shuffled p-value
SETTINGS.PPREGRESSION.whichSpikeInTheta = 'mean';                            % choose between all, first, last and mean
SETTINGS.PPREGRESSION.splitVar = 'all';                                     % choose between {'all'; 'spike'; 'speed'; 'rate'; 'peripheral'; 'tortuosity'}
SETTINGS.PPREGRESSION.hilo = false;                                         % split in 2 according to median 
SETTINGS.PPREGRESSION.PLOT_ON = debugPlot;

%% Construct settings required for resultsSummary
SETTINGS.RESULTSSUMMARY.SHALLOW_THRESHOLD = 90*(pi/180);                    % Define threshold for steep or shallow phase precession
SETTINGS.RESULTSSUMMARY.DIRECTIONALITY_THRESHOLD = 0.1;                     % Define threshold for directionality
SETTINGS.RESULTSSUMMARY.GRIDNESS_THRESHOLD = 0;                             % Define threshold for grid cell