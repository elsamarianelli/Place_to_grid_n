function [ out ] = skaggsInfo( varargin )
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info

if nargin == 0, out = DealWithInputs; return, end;
in = DealWithInputs(varargin{:});

nanmask = ~(isnan(in.rates) | isnan(in.times));
duration = sum(in.times(nanmask));
meanRate = sum(in.rates(nanmask).*in.times(nanmask))/duration;

if meanRate<=0.0
    fprintf(' mean rate is %f \n', meanRate);
    out.bits_per_sec = NaN;
    out.bits_per_spike = NaN;
    return;
end

p_x = in.times ./ duration;
p_r = in.rates ./ meanRate;
% avoid finding 0.log(0) terms in sum
dummy = p_x .* in.rates;
idx = dummy > 0;
out.bitsPerSec = sum(dummy(idx).*log2(p_r(idx)));
out.bitsPerSpike = out.bitsPerSec/meanRate;

end

function in = DealWithInputs(varargin)
defaults.rates = double.empty(64,0);
defaults.times = double.empty(64,0);

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