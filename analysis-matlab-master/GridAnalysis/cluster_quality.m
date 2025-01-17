function [ out ] = cluster_quality ( varargin )
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info

if nargin == 0, out = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

% calculate energy on each electrode
E = sqrt(sum(in.spike(:,:,:).^2,2));

% normalise each waveform by its energy
normdWaves = bsxfun(@rdivide, in.spike, E);

% get principal components - default is 2 prin comps
nSpikes = size(normdWaves,1);
nElectrodes = size(normdWaves,3);
PCA = zeros(nSpikes, in.nComps, nElectrodes);
for i = 1:nElectrodes
    [~, P] = princomp(normdWaves(:,:,i));
    PCA(:, 1:in.nComps, i) = P(:, 1:in.nComps);
end

% get mahalanobis distance
% reshape the PCA matrix
PCA_m = reshape(PCA, nSpikes, prod([in.nComps, nElectrodes]));
M = mahal(PCA_m, PCA_m(in.spikeInd,:));

% get the indices of the spikes not in the cluster
noiseSpkInd = setdiff(1:nSpikes, in.spikeInd);
% calculate L-ratio...
M_noise = M(noiseSpkInd);
% get degrees of freedom
df = prod([in.nComps, nElectrodes]);
% calculate L-ratio
L = sum(1 - chi2cdf(M(noiseSpkInd), df));
L_ratio = L / length(in.spikeInd);

% calculate the isolation distance
if (length(in.spikeInd) < nSpikes/2) % return Nan if more spikes in cluster than 1/2 of total
    [A, ~] = sort(M_noise);
    isolation_dist = A(length(in.spikeInd));
else
    isolation_dist = nan;
end

out.L_ratio = L_ratio;
out.isolation_dist = isolation_dist;
end

function in = DealWithInputs(varargin)
defaults.spike = [];
defaults.spikeInd = [];
defaults.nComps = 2; % number of principal components to calculate
VERSION = 1;

% Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    in = defaults;
    return;
end

if isstruct(varargin{1})
    if nargin > 1 && VERSION ~= varargin{2}
        error(['%s called with version number %g, but code is version %g.\n' ...
            'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
            'the function without a version number.'],mfilename,varargin{2},VERSION);
    end
    in = ModifyExistingFields(defaults,varargin{1});
else
    in = ModifyExistingFields(defaults,varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
